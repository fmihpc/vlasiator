#include <float.h>
#include <getopt.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TAG_LAT 1
#define TAG_BW 2
#define TAG_ACK 3
#define LATBYTES 8
/* BSD sysexits value for a bad command line, so scripts can tell it apart from a failed run. */
#define EXIT_USAGE 64

typedef struct {
    double time;
    int size, window, pingpongs;
} options_t;

/* What one rank measured: latency and forward bandwidth on the link to its
 * right neighbour, and bandwidth back on the link to its left one. */
typedef struct {
    char host[MPI_MAX_PROCESSOR_NAME];
    double lat_sum, lat_min, lat_max;
    long lat_n;
    double fwd_bytes, fwd_time, fwd_worst;
    double back_bytes, back_time, back_worst;
} stats_t;

static int rank, nranks;
static options_t opt;
static char *sbuf, *rbuf;
static MPI_Request *req;

enum { OPT_TIME, OPT_SIZE, OPT_WINDOW, OPT_PINGPONGS, OPT_HELP };

/* In enum order. The parser and the help text both come from this table, so they cannot drift apart. */
/* clang-format off */
static const struct {
    const char *name;
    char letter;
    const char *metavar, *expect, *help;
} option_table[] = {
    {"time", 't', "S", "seconds, more than 0 and at most 604800", "seconds to measure after one warm-up round (default 60)"},
    {"size", 's', "B", "bytes from 1 to 1073741824", "bytes per bandwidth message (default 1048576)"},
    {"window", 'w', "N", "a whole number from 1 to 1024", "bandwidth messages in flight at once (default 16; window x size at most 1 GiB)"},
    {"pingpongs", 'p', "N", "a whole number from 1 to 1000000", "latency round trips per link per round (default 100)"},
    {"help", 'h', NULL, NULL, "show this help and exit"},
};
/* clang-format on */
#define NOPT ((int)(sizeof option_table / sizeof option_table[0]))

static void usage(FILE *f, const char *prog) {
    int k;

    fprintf(f, "usage: %s [options]\n", prog);
    fprintf(f, "Each rank exchanges messages with its ring neighbours for a set time, then rank 0 reports latency and bandwidth per link and on average.\n");
    for (k = 0; k < NOPT; ++k) {
        const char *mv = option_table[k].metavar;
        int pad = 14 - (int)(strlen(option_table[k].name) + (mv ? 1 + strlen(mv) : 0));

        fprintf(f, "  -%c, --%s%s%s%*s %s\n", option_table[k].letter, option_table[k].name, mv ? "=" : "", mv ? mv : "", pad, "", option_table[k].help);
    }
    fprintf(f, "Latency is half the round trip of %d-byte messages. Run one rank per node to measure only links between nodes.\n", LATBYTES);
    fprintf(f, "Exit status: 0 measured, %d bad command line.\n", EXIT_USAGE);
}

static int usage_hint(const char *prog) {
    if (rank == 0) {
        fprintf(stderr, "Try '%s --help'.\n", prog);
    }
    return EXIT_USAGE;
}

/* The whole string must be a number in (lo, hi]. NaN fails every comparison, inf fails the upper bound. */
static int parse_double(const char *s, double lo, double hi, double *v) {
    char *end;
    double x = strtod(s, &end);

    if (end == s || *end != '\0' || !(x > lo) || !(x <= hi)) {
        return 0;
    }
    *v = x;
    return 1;
}

static int parse_int(const char *s, long lo, long hi, int *v) {
    char *end;
    long x = strtol(s, &end, 10);

    /* strtol saturates on overflow, which the range check then rejects. */
    if (end == s || *end != '\0' || x < lo || x > hi) {
        return 0;
    }
    *v = (int)x;
    return 1;
}

static int option_index(int letter) {
    int k;

    for (k = 0; k < NOPT; ++k) {
        if (option_table[k].letter == letter) {
            return k;
        }
    }
    return -1;
}

/* Returns 0 to go on, -1 after --help, or EXIT_USAGE. Every rank parses the same arguments, so only rank 0 reports. */
static int parse_args(int argc, char **argv, options_t *o) {
    const char *prog = argc > 0 ? argv[0] : "ringperf";
    struct option longopts[NOPT + 1];
    char shortopts[2 * NOPT + 1];
    int c, k, n = 0, ok;

    o->time = 60.0;
    o->size = 1048576;
    o->window = 16;
    o->pingpongs = 100;
    for (k = 0; k < NOPT; ++k) {
        longopts[k].name = option_table[k].name;
        longopts[k].has_arg = option_table[k].metavar ? required_argument : no_argument;
        longopts[k].flag = NULL;
        longopts[k].val = option_table[k].letter;
        shortopts[n++] = option_table[k].letter;
        if (option_table[k].metavar) {
            shortopts[n++] = ':';
        }
    }
    memset(&longopts[NOPT], 0, sizeof longopts[NOPT]);
    shortopts[n] = '\0';
    opterr = rank == 0;
    while ((c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
        k = option_index(c);
        if (k < 0) {
            /* getopt_long has already said what was wrong. */
            return usage_hint(prog);
        }
        ok = 1;
        switch (k) {
        case OPT_TIME:
            ok = parse_double(optarg, 0.0, 604800.0, &o->time);
            break;
        case OPT_SIZE:
            ok = parse_int(optarg, 1, 1073741824, &o->size);
            break;
        case OPT_WINDOW:
            ok = parse_int(optarg, 1, 1024, &o->window);
            break;
        case OPT_PINGPONGS:
            ok = parse_int(optarg, 1, 1000000, &o->pingpongs);
            break;
        case OPT_HELP:
            if (rank == 0) {
                usage(stdout, prog);
            }
            return -1;
        }
        if (!ok) {
            if (rank == 0) {
                fprintf(stderr, "%s: invalid value '%s' for --%s: expected %s\n", prog, optarg, option_table[k].name, option_table[k].expect);
            }
            return usage_hint(prog);
        }
    }
    if (optind < argc) {
        if (rank == 0) {
            fprintf(stderr, "%s: unexpected argument '%s'\n", prog, argv[optind]);
        }
        return usage_hint(prog);
    }
    /* Every message in a window needs its own receive buffer. */
    if ((double)o->window * o->size > 1073741824.0) {
        if (rank == 0) {
            fprintf(stderr, "%s: --window x --size is %.0f bytes, more than 1 GiB\n", prog, (double)o->window * o->size);
        }
        return usage_hint(prog);
    }
    return 0;
}

/* Two ranks share a single link; a lone rank has none. */
static int nlinks(void) {
    return nranks >= 3 ? nranks : nranks - 1;
}

/* Alternating groups keep every rank on at most one link at a time, so no
 * rank's timings include work on another link. An odd ring needs a third group
 * for the link that closes it. */
static int link_group(int k) {
    return nranks % 2 == 1 && k == nranks - 1 ? 2 : k % 2;
}

static int ngroups(void) {
    return nranks >= 3 ? 2 + nranks % 2 : nlinks();
}

static void lat_initiate(int peer, stats_t *st, int record) {
    char b[LATBYTES];
    double t;
    int i;

    memset(b, 0, sizeof b);
    /* Untimed first exchange: the partner may still be busy with its previous link. */
    MPI_Send(b, LATBYTES, MPI_BYTE, peer, TAG_LAT, MPI_COMM_WORLD);
    MPI_Recv(b, LATBYTES, MPI_BYTE, peer, TAG_LAT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (i = 0; i < opt.pingpongs; ++i) {
        t = MPI_Wtime();
        MPI_Send(b, LATBYTES, MPI_BYTE, peer, TAG_LAT, MPI_COMM_WORLD);
        MPI_Recv(b, LATBYTES, MPI_BYTE, peer, TAG_LAT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /* Half the round trip: a one-way time would need synchronized clocks. */
        t = (MPI_Wtime() - t) / 2;
        if (record) {
            st->lat_sum += t;
            ++st->lat_n;
            if (t < st->lat_min) {
                st->lat_min = t;
            }
            if (t > st->lat_max) {
                st->lat_max = t;
            }
        }
    }
}

static void lat_answer(int peer) {
    char b[LATBYTES];
    int i;

    for (i = 0; i <= opt.pingpongs; ++i) {
        MPI_Recv(b, LATBYTES, MPI_BYTE, peer, TAG_LAT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(b, LATBYTES, MPI_BYTE, peer, TAG_LAT, MPI_COMM_WORLD);
    }
}

static void bw_send(int peer, double *bytes, double *time, double *worst, int record) {
    double t = MPI_Wtime(), sent = (double)opt.window * opt.size;
    int i, ack;

    for (i = 0; i < opt.window; ++i) {
        MPI_Isend(sbuf, opt.size, MPI_BYTE, peer, TAG_BW, MPI_COMM_WORLD, &req[i]);
    }
    MPI_Waitall(opt.window, req, MPI_STATUSES_IGNORE);
    /* Local completion only means the data left the buffer; the acknowledgement means it arrived. */
    MPI_Recv(&ack, 1, MPI_INT, peer, TAG_ACK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    t = MPI_Wtime() - t;
    if (record && t > 0) {
        *bytes += sent;
        *time += t;
        if (sent / t < *worst) {
            *worst = sent / t;
        }
    }
}

static void bw_recv(int peer) {
    int i, ack = 0;

    /* One buffer per message: concurrent receives into the same buffer would be erroneous. */
    for (i = 0; i < opt.window; ++i) {
        MPI_Irecv(rbuf + (size_t)i * opt.size, opt.size, MPI_BYTE, peer, TAG_BW, MPI_COMM_WORLD, &req[i]);
    }
    MPI_Waitall(opt.window, req, MPI_STATUSES_IGNORE);
    MPI_Send(&ack, 1, MPI_INT, peer, TAG_ACK, MPI_COMM_WORLD);
}

static void one_round(stats_t *st, int record) {
    int own = rank < nlinks() ? rank : -1, left = (rank - 1 + nranks) % nranks, g;

    for (g = 0; g < ngroups(); ++g) {
        if (own >= 0 && link_group(own) == g) {
            /* The left end of a link times latency and the forward direction, the right end the way back. */
            lat_initiate((rank + 1) % nranks, st, record);
            bw_send((rank + 1) % nranks, &st->fwd_bytes, &st->fwd_time, &st->fwd_worst, record);
            bw_recv((rank + 1) % nranks);
        } else if (left < nlinks() && left != own && link_group(left) == g) {
            lat_answer(left);
            bw_recv(left);
            bw_send(left, &st->back_bytes, &st->back_time, &st->back_worst, record);
        }
    }
}

static void report(const stats_t *all, int rounds, double elapsed) {
    double lat[2] = {0.0, 0.0}, bw[2] = {0.0, 0.0}, low[2] = {DBL_MAX, DBL_MAX};
    int cnt[2] = {0, 0}, lowlink[2] = {0, 0}, k, c;

    printf("ringperf: %d timed round(s) in %.1f s after one warm-up round\n", rounds, elapsed);
    for (k = 0; k < nlinks(); ++k) {
        const stats_t *a = &all[k], *b = &all[(k + 1) % nranks];
        double l = a->lat_n > 0 ? a->lat_sum / a->lat_n : 0.0, f = a->fwd_time > 0 ? a->fwd_bytes / a->fwd_time : 0.0, r = b->back_time > 0 ? b->back_bytes / b->back_time : 0.0;

        c = strcmp(a->host, b->host) != 0;
        printf("  link %d: rank %d (%s) - rank %d (%s): latency %.2f us (min %.2f, max %.2f); bandwidth %.3f GB/s forward (worst round %.3f), %.3f GB/s back (worst round %.3f)\n", k, k, a->host, (k + 1) % nranks, b->host, 1e6 * l, 1e6 * a->lat_min, 1e6 * a->lat_max, 1e-9 * f, 1e-9 * a->fwd_worst, 1e-9 * r, 1e-9 * b->back_worst);
        ++cnt[c];
        lat[c] += l;
        bw[c] += f + r;
        if ((f < r ? f : r) < low[c]) {
            low[c] = f < r ? f : r;
            lowlink[c] = k;
        }
    }
    for (c = 1; c >= 0; --c) {
        if (cnt[c] > 0) {
            printf("ringperf: %s: %d link(s), average latency %.2f us, average bandwidth %.3f GB/s, lowest %.3f GB/s on link %d\n", c ? "between hosts" : "within hosts", cnt[c], 1e6 * lat[c] / cnt[c], 1e-9 * bw[c] / (2 * cnt[c]), 1e-9 * low[c], lowlink[c]);
        }
    }
}

int main(int argc, char **argv) {
    stats_t st, *all = NULL;
    double t0, elapsed = 0.0;
    int len, rounds = 0, stop = 0, code;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    code = parse_args(argc, argv, &opt);
    if (code != 0) {
        MPI_Finalize();
        return code < 0 ? 0 : EXIT_USAGE;
    }
    if (nranks < 2) {
        if (rank == 0) {
            printf("ringperf: a single rank has no links to measure\n");
        }
        MPI_Finalize();
        return 0;
    }

    memset(&st, 0, sizeof st);
    st.lat_min = st.fwd_worst = st.back_worst = DBL_MAX;
    MPI_Get_processor_name(st.host, &len);
    sbuf = (char *)malloc((size_t)opt.size);
    rbuf = (char *)malloc((size_t)opt.window * opt.size);
    req = (MPI_Request *)malloc((size_t)opt.window * sizeof *req);
    if (rank == 0) {
        all = (stats_t *)malloc((size_t)nranks * sizeof *all);
    }
    if (!sbuf || !rbuf || !req || (rank == 0 && !all)) {
        fprintf(stderr, "ringperf: rank %d: out of memory\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    /* Touch every page now, so first-touch page faults don't land in the timings. */
    memset(sbuf, 1, (size_t)opt.size);
    memset(rbuf, 0, (size_t)opt.window * opt.size);
    if (rank == 0) {
        printf("ringperf: %d ranks, %d link(s); latency from %d round trips of %d bytes, bandwidth from %d x %d-byte messages in flight, per link and round, both directions; %.0f s\n", nranks, nlinks(), opt.pingpongs, LATBYTES, opt.window, opt.size, opt.time);
        fflush(stdout);
    }

    /* Untimed, because first contact sets up connections and memory registration, which would skew the first samples. */
    one_round(&st, 0);
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    while (!stop) {
        one_round(&st, 1);
        ++rounds;
        if (rank == 0) {
            elapsed = MPI_Wtime() - t0;
            stop = elapsed >= opt.time;
        }
        /* Rank 0's clock decides for everyone, so no rank starts a round its partners won't join. */
        MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    MPI_Gather(&st, (int)sizeof st, MPI_BYTE, all, (int)sizeof st, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        report(all, rounds, elapsed);
    }
    MPI_Finalize();
    return 0;
}
