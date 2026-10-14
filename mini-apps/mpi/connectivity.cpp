#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TAG_PROBE 1
#define TAG_REPORT 2
#define TAG_GO 3
#define TAG_LAT 4
#define TAG_HDR 5
#define TAG_STATS 6
#define MAXRANGES 32
#define TOPN 10
#define TOPPAIRS 3
#define MAXEV 256
#define TAG_EV 7

enum { PENDING = 0, DONE = 1, FAILED = 2 };

typedef struct {
    char host[MPI_MAX_PROCESSOR_NAME];
    int nfail;
    int nranges;
    int range[MAXRANGES][2];
} report_t;

/* Round trips of pair (a, b), timed and reported by a < b. */
typedef struct {
    int a, b, n, nspike, stalled;
    float min, med, max;
} pairstat_t;

/* One spike. Steps run in lockstep, so spikes sharing round and step happened at the same moment. */
typedef struct {
    int round, step, a, b;
    float rtt;
} event_t;

typedef struct {
    long hp, n;
} hpcount_t;

static int rank, size;
static double timeout;
static char host[MPI_MAX_PROCESSOR_NAME];

/* In-flight requests, each with the status byte it settles. */
static MPI_Request *req;
static char **flag;
static int *idx;
static MPI_Status *st;
static int nact;

static void post(int send, void *buf, int count, MPI_Datatype type, int peer, int tag, char *f) {
    int err = send ? MPI_Isend(buf, count, type, peer, tag, MPI_COMM_WORLD, &req[nact])
                   : MPI_Irecv(buf, count, type, peer, tag, MPI_COMM_WORLD, &req[nact]);
    if (err != MPI_SUCCESS) {
        *f = FAILED;
        return;
    }
    flag[nact++] = f;
}

static void progress(void) {
    int k, n = 0, out = 0, err;

    if (nact == 0) {
        return;
    }
    err = MPI_Testsome(nact, req, &out, idx, st);
    if (err != MPI_SUCCESS && err != MPI_ERR_IN_STATUS) {
        /* Outcome per request is unknown; treat all as failed rather than spin. */
        for (k = 0; k < nact; ++k) {
            *flag[k] = FAILED;
        }
        nact = 0;
        return;
    }
    if (out == MPI_UNDEFINED) {
        out = 0;
    }
    for (k = 0; k < out; ++k) {
        *flag[idx[k]] = (err == MPI_ERR_IN_STATUS && st[k].MPI_ERROR != MPI_SUCCESS) ? FAILED : DONE;
    }
    for (k = 0; k < nact; ++k) {
        if (*flag[k] == PENDING) {
            req[n] = req[k];
            flag[n] = flag[k];
            ++n;
        }
    }
    nact = n;
}

/* Polls until all n requests complete (1), or the deadline passes or an error occurs (0). */
static int wait_until(MPI_Request *r, int n, double deadline) {
    int done;

    for (;;) {
        if (MPI_Testall(n, r, &done, MPI_STATUSES_IGNORE) != MPI_SUCCESS) {
            return 0;
        }
        if (done) {
            return 1;
        }
        if (MPI_Wtime() > deadline) {
            return 0;
        }
    }
}

static int next_run(const char *s, int n, int *j, int *lo, int *hi) {
    while (*j < n && s[*j] == DONE) {
        ++*j;
    }
    if (*j == n) {
        return 0;
    }
    *lo = *j;
    while (*j < n && s[*j] != DONE) {
        ++*j;
    }
    *hi = *j - 1;
    return 1;
}

static void print_runs(const char *s, int n) {
    int j = 0, lo, hi, first = 1;

    while (next_run(s, n, &j, &lo, &hi)) {
        printf("%s%d", first ? "" : ",", lo);
        if (hi > lo) {
            printf("-%d", hi);
        }
        first = 0;
    }
}

static int count_not_done(const char *s, int n) {
    int j, c = 0;

    for (j = 0; j < n; ++j) {
        c += s[j] != DONE;
    }
    return c;
}

static report_t *g_rep;
static const double *g_key;

static int by_host(const void *a, const void *b) {
    return strcmp(g_rep[*(const int *)a].host, g_rep[*(const int *)b].host);
}

static int by_key_desc(const void *a, const void *b) {
    double x = g_key[*(const long *)a], y = g_key[*(const long *)b];
    return (x < y) - (x > y);
}

static int by_float(const void *a, const void *b) {
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

static const int *g_hid;
static int g_nh;
static const pairstat_t *g_ps;

/* Unordered host pair as one number, so rank pairs and spikes can be grouped by the hosts they connect. */
static long hostpair(int a, int b) {
    int x = g_hid[a], y = g_hid[b];

    return x < y ? (long)x * g_nh + y : (long)y * g_nh + x;
}

static int by_hostpair_med(const void *a, const void *b) {
    const pairstat_t *x = &g_ps[*(const long *)a], *y = &g_ps[*(const long *)b];
    long hx = hostpair(x->a, x->b), hy = hostpair(y->a, y->b);

    if (hx != hy) {
        return (hx > hy) - (hx < hy);
    }
    return (x->med > y->med) - (x->med < y->med);
}

static int by_max_desc(const void *a, const void *b) {
    float x = g_ps[*(const long *)a].max, y = g_ps[*(const long *)b].max;

    return (x < y) - (x > y);
}

static int by_rtt_desc(const void *a, const void *b) {
    float x = ((const event_t *)a)->rtt, y = ((const event_t *)b)->rtt;

    return (x < y) - (x > y);
}

static int by_moment(const void *a, const void *b) {
    const event_t *x = (const event_t *)a, *y = (const event_t *)b;
    long hx, hy;

    if (x->round != y->round) {
        return (x->round > y->round) - (x->round < y->round);
    }
    if (x->step != y->step) {
        return (x->step > y->step) - (x->step < y->step);
    }
    hx = hostpair(x->a, x->b);
    hy = hostpair(y->a, y->b);
    return (hx > hy) - (hx < hy);
}

static int by_count_desc(const void *a, const void *b) {
    long x = ((const hpcount_t *)a)->n, y = ((const hpcount_t *)b)->n;

    return (x < y) - (x > y);
}

static int by_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* Circle method: with P even, every pair meets exactly once in steps 0..P-2; a partner >= size is a bye. */
static int partner(int i, int step, int P) {
    if (i == P - 1) {
        return step;
    }
    if (i == step) {
        return P - 1;
    }
    return ((2 * step - i) % (P - 1) + (P - 1)) % (P - 1);
}

static int lat_recv(int *buf, int p, MPI_Request *r) {
    return MPI_Irecv(buf, 1, MPI_INT, p, TAG_LAT, MPI_COMM_WORLD, r) == MPI_SUCCESS;
}

static int lat_send(int p, MPI_Request *r) {
    static int one = 1;
    return MPI_Isend(&one, 1, MPI_INT, p, TAG_LAT, MPI_COMM_WORLD, r) == MPI_SUCCESS;
}

/* One timed round trip per pair per round. Returns the number of stalls; sets
 * *gave_up if a barrier timed out. A stalled pair is dropped for good, because
 * its abandoned requests would match later messages. */
static int latency(int rounds, float *samp, int *pstep, char *dead, int *inbuf, int *gave_up) {
    int P = size + (size & 1), r, step, stalls = 0;
    MPI_Request rq[2];

    *gave_up = 0;
    for (r = 0; r < rounds; ++r) {
        for (step = 0; step < P - 1; ++step) {
            int p = partner(rank, step, P), ok;
            double dl, t;

            /* Lockstep caps the delay a stalled pair imposes on everyone else at one timeout, so stalls can't cascade. */
            if (MPI_Ibarrier(MPI_COMM_WORLD, &rq[0]) != MPI_SUCCESS || !wait_until(rq, 1, MPI_Wtime() + 2 * timeout)) {
                printf("rank %d (%s): latency barrier timed out in round %d, giving up\n", rank, host, r);
                fflush(stdout);
                *gave_up = 1;
                return stalls;
            }
            if (p >= size || dead[p]) {
                continue;
            }
            dl = MPI_Wtime() + timeout;
            if (rank < p) {
                /* Untimed first exchange absorbs the partner leaving the barrier later than us. */
                ok = lat_recv(&inbuf[p], p, &rq[0]) &&
                     lat_send(p, &rq[1]) &&
                     wait_until(rq, 2, dl) &&
                     lat_recv(&inbuf[p], p, &rq[0]);
                if (ok) {
                    t = MPI_Wtime();
                    ok = lat_send(p, &rq[1]) && wait_until(rq, 2, dl);
                    if (ok) {
                        int i = p - rank - 1;

                        /* Indexed by round, and each pair keeps its step every round, which places every spike in time. */
                        samp[(size_t)i * rounds + r] = (float)(MPI_Wtime() - t);
                        pstep[i] = step;
                    }
                }
            } else {
                /* Preposting the timed ping's receive keeps it off the unexpected-message path. */
                ok = lat_recv(&inbuf[p], p, &rq[0]) &&
                     wait_until(rq, 1, dl) &&
                     lat_recv(&inbuf[p], p, &rq[0]) &&
                     lat_send(p, &rq[1]) &&
                     wait_until(&rq[1], 1, dl) &&
                     wait_until(rq, 1, dl) &&
                     lat_send(p, &rq[1]) &&
                     wait_until(&rq[1], 1, dl);
            }
            if (!ok) {
                dead[p] = 1;
                ++stalls;
                printf("rank %d (%s): latency exchange with rank %d stalled in round %d\n", rank, host, p, r);
                fflush(stdout);
            }
        }
    }
    return stalls;
}

static double spike_threshold(double med, double factor) {
    double tick = MPI_Wtick();

    /* Clamp to timer resolution so a zero median doesn't make every sample a spike. */
    return factor * (med > tick ? med : tick);
}

/* Fills out[] and copies this rank's largest spikes, at most MAXEV, into ev; returns how many. Missing samples are negative. */
static int pair_stats(int rounds, double factor, const float *samp, const int *pstep, const char *dead, pairstat_t *out, int nmine, event_t *ev) {
    float *x = (float *)malloc(((size_t)rounds + 1) * sizeof *x);
    event_t *all;
    long nall = 0, k = 0;
    int i, r, n;

    if (!x) {
        fprintf(stderr, "rank %d: out of memory for latency statistics\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (i = 0; i < nmine; ++i) {
        pairstat_t *s = &out[i];
        const float *y = samp + (size_t)i * rounds;
        double thr;

        memset(s, 0, sizeof *s);
        s->a = rank;
        s->b = rank + 1 + i;
        s->stalled = dead[s->b];
        for (n = 0, r = 0; r < rounds; ++r) {
            if (y[r] >= 0) {
                x[n++] = y[r];
            }
        }
        s->n = n;
        if (n == 0) {
            continue;
        }
        qsort(x, n, sizeof *x, by_float);
        s->min = x[0];
        s->med = x[n / 2];
        s->max = x[n - 1];
        thr = spike_threshold(s->med, factor);
        for (r = 0; r < n; ++r) {
            s->nspike += x[r] > thr;
        }
        nall += s->nspike;
    }
    all = (event_t *)malloc(((size_t)nall + 1) * sizeof *all);
    if (!all) {
        fprintf(stderr, "rank %d: out of memory for latency statistics\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (i = 0; i < nmine; ++i) {
        const float *y = samp + (size_t)i * rounds;
        double thr;

        if (out[i].n == 0) {
            continue;
        }
        thr = spike_threshold(out[i].med, factor);
        for (r = 0; r < rounds; ++r) {
            if (y[r] >= 0 && y[r] > thr) {
                all[k].round = r;
                all[k].step = pstep[i];
                all[k].a = rank;
                all[k].b = out[i].b;
                all[k].rtt = y[r];
                ++k;
            }
        }
    }
    qsort(all, k, sizeof *all, by_rtt_desc);
    n = k < MAXEV ? (int)k : MAXEV;
    memcpy(ev, all, (size_t)n * sizeof *ev);
    free(x);
    free(all);
    return n;
}

static void print_pair(const pairstat_t *s) {
    printf("    %d (%s) - %d (%s): ", s->a, g_rep[s->a].host, s->b, g_rep[s->b].host);
}

/* Prints the latency report on rank 0; returns 1 if any outlier was flagged. */
static int analyze(const pairstat_t *ps, long nps, event_t *ev, long nev, int ntrunc, int rounds, double factor, double elapsed) {
    int *order = (int *)malloc((size_t)size * sizeof *order), *hid = (int *)malloc((size_t)size * sizeof *hid);
    long *sel = (long *)malloc(((size_t)nps + size) * sizeof *sel), *hsp, *hsm, nm[2], nsel, nslow, p;
    long nsamples = 0, nspikes = 0;
    float *meds = (float *)malloc(((size_t)nps + 1) * sizeof *meds);
    double *key = (double *)malloc(((size_t)nps + size) * sizeof *key), *rate;
    double typ[2], thr, tick = MPI_Wtick();
    const char **hname = (const char **)malloc((size_t)size * sizeof *hname);
    int nh = 0, j, k, c;
    long q, e, nrun, *rstart, *rlen, *rorder;

    if (!order || !hid || !sel || !meds || !key || !hname) {
        fprintf(stderr, "rank 0: out of memory for latency statistics\n");
        return 0;
    }
    for (j = 0; j < size; ++j) {
        order[j] = j;
    }
    qsort(order, size, sizeof *order, by_host);
    for (k = 0; k < size; ++k) {
        if (k == 0 || strcmp(g_rep[order[k]].host, g_rep[order[k - 1]].host) != 0) {
            hname[nh++] = g_rep[order[k]].host;
        }
        hid[order[k]] = nh - 1;
    }
    g_hid = hid;
    g_nh = nh;

    /* Same-host and inter-host round trips differ by an order of magnitude, so each gets its own baseline. */
    for (c = 0; c < 2; ++c) {
        nm[c] = 0;
        for (p = 0; p < nps; ++p) {
            if (ps[p].n > 0 && (hid[ps[p].a] != hid[ps[p].b]) == c) {
                meds[nm[c]++] = ps[p].med;
            }
        }
        qsort(meds, nm[c], sizeof *meds, by_float);
        /* Lower median: with 4 hosts, exactly half the inter-host pairs touch any one host. */
        typ[c] = nm[c] > 0 ? meds[(nm[c] - 1) / 2] : 0.0;
        if (typ[c] < tick) {
            typ[c] = tick;
        }
    }
    for (p = 0; p < nps; ++p) {
        nsamples += ps[p].n;
        nspikes += ps[p].nspike;
    }
    printf("latency: %d round(s) of %d step(s), %ld timed round trips in %.1f s, timer resolution %.3g us\n", rounds, size + (size & 1) - 1, nsamples, elapsed, 1e6 * tick);
    for (c = 0; c < 2; ++c) {
        if (nm[c] > 0) {
            printf("  %s: %ld pair(s), typical round trip %.2f us (median of pair medians)\n", c ? "different hosts" : "same host", nm[c], 1e6 * typ[c]);
        }
    }

    nsel = 0;
    for (p = 0; p < nps; ++p) {
        if (ps[p].n == 0) {
            continue;
        }
        key[p] = ps[p].med / typ[hid[ps[p].a] != hid[ps[p].b]];
        if (key[p] > factor) {
            sel[nsel++] = p;
        }
    }
    g_key = key;
    qsort(sel, nsel, sizeof *sel, by_key_desc);
    nslow = nsel;
    printf("  %ld slow pair(s): median > %gx the typical round trip%s\n", nslow, factor, nslow ? "; worst:" : "");
    for (k = 0; k < nsel && k < TOPN; ++k) {
        print_pair(&ps[sel[k]]);
        printf("median %.2f us (%.0fx), min %.2f us, max %.2f us\n", 1e6 * ps[sel[k]].med, key[sel[k]], 1e6 * ps[sel[k]].min, 1e6 * ps[sel[k]].max);
    }

    nsel = 0;
    for (p = 0; p < nps; ++p) {
        nsel += ps[p].nspike > 0;
    }
    printf("  %ld spike(s) (round trip > %gx its pair's median) in %ld pair(s)\n", nspikes, factor, nsel);

    /* NICs, cables and switch ports belong to hosts, and which ranks happened to be paired is incidental, so report by host pair. */
    nsel = 0;
    for (p = 0; p < nps; ++p) {
        if (ps[p].n > 0) {
            sel[nsel++] = p;
        }
    }
    rstart = (long *)malloc(((size_t)nsel + 1) * sizeof *rstart);
    rlen = (long *)malloc(((size_t)nsel + 1) * sizeof *rlen);
    rorder = (long *)malloc(((size_t)nsel + 1) * sizeof *rorder);
    if (!rstart || !rlen || !rorder) {
        fprintf(stderr, "rank 0: out of memory for host pair statistics\n");
        return nslow > 0;
    }
    g_ps = ps;
    qsort(sel, nsel, sizeof *sel, by_hostpair_med);
    for (nrun = 0, q = 0; q < nsel; q = e) {
        long hp = hostpair(ps[sel[q]].a, ps[sel[q]].b);
        double mx = 0.0;

        for (e = q; e < nsel && hostpair(ps[sel[e]].a, ps[sel[e]].b) == hp; ++e) {
            if (ps[sel[e]].max > mx) {
                mx = ps[sel[e]].max;
            }
        }
        rstart[nrun] = q;
        rlen[nrun] = e - q;
        key[nrun] = mx;
        rorder[nrun] = nrun;
        ++nrun;
    }
    g_key = key;
    qsort(rorder, nrun, sizeof *rorder, by_key_desc);
    printf("  %ld host pair(s), largest round trip first", nrun);
    if (nrun > TOPN) {
        printf(" (top %d)", TOPN);
    }
    printf("\n");
    for (k = 0; k < nrun && k < TOPN; ++k) {
        long r0 = rstart[rorder[k]], rn = rlen[rorder[k]], hp = hostpair(ps[sel[r0]].a, ps[sel[r0]].b), nrt = 0, nsp = 0;
        /* Each run is sorted by median, so this is its lower median. */
        double typical = ps[sel[r0 + (rn - 1) / 2]].med;

        for (q = r0; q < r0 + rn; ++q) {
            nrt += ps[sel[q]].n;
            nsp += ps[sel[q]].nspike;
        }
        printf("    %s - %s: %ld rank pair(s), %ld round trips, %ld spike(s), typical %.2f us, max %.2f us\n", hname[hp / nh], hname[hp % nh], rn, nrt, nsp, 1e6 * typical, 1e6 * key[rorder[k]]);
        if (nsp == 0) {
            continue;
        }
        /* Still show rank pairs: ranks bound to different NICs of one node can behave differently, which host totals would hide. */
        qsort(sel + r0, rn, sizeof *sel, by_max_desc);
        for (q = r0; q < r0 + rn && q < r0 + TOPPAIRS; ++q) {
            printf("      %d - %d: max %.2f us, median %.2f us, %d spike(s)\n", ps[sel[q]].a, ps[sel[q]].b, 1e6 * ps[sel[q]].max, 1e6 * ps[sel[q]].med, ps[sel[q]].nspike);
        }
    }

    if (nev > 0) {
        long ng = 0, *gs = (long *)malloc((size_t)nev * sizeof *gs), *gl = (long *)malloc((size_t)nev * sizeof *gl), *go = (long *)malloc((size_t)nev * sizeof *go);
        double *gk = (double *)malloc((size_t)nev * sizeof *gk);
        hpcount_t *hc = (hpcount_t *)malloc((size_t)nev * sizeof *hc);

        if (!gs || !gl || !go || !gk || !hc) {
            fprintf(stderr, "rank 0: out of memory for spike events\n");
            return nslow > 0;
        }
        qsort(ev, nev, sizeof *ev, by_moment);
        for (q = 0; q < nev; q = e) {
            double mx = 0.0;

            for (e = q; e < nev && ev[e].round == ev[q].round && ev[e].step == ev[q].step; ++e) {
                if (ev[e].rtt > mx) {
                    mx = ev[e].rtt;
                }
            }
            gs[ng] = q;
            gl[ng] = e - q;
            gk[ng] = mx;
            go[ng] = ng;
            ++ng;
        }
        g_key = gk;
        qsort(go, ng, sizeof *go, by_key_desc);
        printf("  spikes by moment (same round and step), largest first");
        if (ng > TOPN) {
            printf(" (top %d)", TOPN);
        }
        printf("\n");
        for (k = 0; k < ng && k < TOPN; ++k) {
            long g0 = gs[go[k]], gn = gl[go[k]], nhc = 0;
            double mn = ev[g0].rtt;

            /* Within a moment, events are sorted by host pair, so equal host pairs are adjacent. */
            for (q = g0; q < g0 + gn; ++q) {
                long hp = hostpair(ev[q].a, ev[q].b);

                if (ev[q].rtt < mn) {
                    mn = ev[q].rtt;
                }
                if (nhc == 0 || hc[nhc - 1].hp != hp) {
                    hc[nhc].hp = hp;
                    hc[nhc].n = 0;
                    ++nhc;
                }
                ++hc[nhc - 1].n;
            }
            qsort(hc, nhc, sizeof *hc, by_count_desc);
            if (gn == 1) {
                printf("    round %d, step %d: 1 spike, %.2f us:", ev[g0].round, ev[g0].step, 1e6 * mn);
            } else {
                printf("    round %d, step %d: %ld spikes, %.2f to %.2f us:", ev[g0].round, ev[g0].step, gn, 1e6 * mn, 1e6 * gk[go[k]]);
            }
            for (q = 0; q < nhc && q < TOPPAIRS; ++q) {
                printf("%s %s - %s (%ld)", q ? "," : "", hname[hc[q].hp / nh], hname[hc[q].hp % nh], hc[q].n);
            }
            if (nhc > TOPPAIRS) {
                printf(", +%ld more host pair(s)", nhc - TOPPAIRS);
            }
            printf("\n");
        }
        if (ntrunc > 0) {
            printf("    each rank reports at most its %d largest spikes and %d rank(s) had more, so these counts are lower bounds\n", MAXEV, ntrunc);
        }
    }

    /* A host is flagged when most of its links to other hosts spike: a faulty
     * NIC raises all of them, while a host merely talking to a faulty one has
     * a minority of bad links. The lower median makes "most" a strict
     * majority, so with 3 hosts an innocent host with one bad partner is not
     * flagged. */
    nsel = 0;
    if (nh >= 3) {
        long nl = 0;
        double base, mtyp;

        hsp = (long *)calloc((size_t)nh * nh, sizeof *hsp);
        hsm = (long *)calloc((size_t)nh * nh, sizeof *hsm);
        rate = (double *)malloc(((size_t)nh * nh + 1) * sizeof *rate);
        if (!hsp || !hsm || !rate) {
            fprintf(stderr, "rank 0: out of memory for host statistics\n");
            return nslow > 0;
        }
        for (p = 0; p < nps; ++p) {
            int ha = hid[ps[p].a], hb = hid[ps[p].b];

            if (ps[p].n == 0 || ha == hb) {
                continue;
            }
            hsp[(size_t)ha * nh + hb] += ps[p].nspike;
            hsm[(size_t)ha * nh + hb] += ps[p].n;
            hsp[(size_t)hb * nh + ha] += ps[p].nspike;
            hsm[(size_t)hb * nh + ha] += ps[p].n;
        }
        for (j = 0; j < nh; ++j) {
            for (k = j + 1; k < nh; ++k) {
                if (hsm[(size_t)j * nh + k] > 0) {
                    rate[nl++] = (double)hsm[(size_t)j * nh + k];
                }
            }
        }
        qsort(rate, nl, sizeof *rate, by_double);
        mtyp = nl > 0 ? rate[nl / 2] : 1.0;
        for (j = 0; j < nh; ++j) {
            nl = 0;
            for (k = 0; k < nh; ++k) {
                if (k != j && hsm[(size_t)j * nh + k] > 0) {
                    rate[nl++] = (double)hsp[(size_t)j * nh + k] / hsm[(size_t)j * nh + k];
                }
            }
            qsort(rate, nl, sizeof *rate, by_double);
            key[j] = nl > 0 ? rate[(nl - 1) / 2] : -1.0;
        }
        /* Baseline from host scores, not from links: with few hosts, most links can touch one bad host. */
        nl = 0;
        for (j = 0; j < nh; ++j) {
            if (key[j] >= 0) {
                rate[nl++] = key[j];
            }
        }
        qsort(rate, nl, sizeof *rate, by_double);
        base = nl > 0 ? rate[(nl - 1) / 2] : 0.0;
        /* One extra spike per typical link keeps lone spikes on a quiet machine from flagging anything. */
        thr = factor * (base + 1.0 / mtyp);
        for (j = 0; j < nh; ++j) {
            if (key[j] > thr) {
                sel[nsel++] = j;
            }
        }
        g_key = key;
        qsort(sel, nsel, sizeof *sel, by_key_desc);
        printf("  %ld host(s) where most links to other hosts spike > %gx as often as on a typical host%s\n", nsel, factor, nsel ? ":" : "");
        for (k = 0; k < nsel; ++k) {
            printf("    %s: most links spike at a rate >= %.3g; typical host %.3g\n", hname[sel[k]], key[sel[k]], base);
        }
    }
    fflush(stdout);
    return nslow > 0 || nsel > 0;
}

int main(int argc, char **argv) {
    int len, j, k, lo, hi, nfail, nsend, rounds = 0, nmine, hdr[3], bad;
    double t0, w, factor = 10.0, tl, tl_end;
    size_t z;
    char *rstate, *sstate, *end, rep = PENDING, *dead;
    int *rbuf, *pstep, *inbuf;
    float *samp;
    pairstat_t *mystats;
    event_t *ev;
    report_t mine;

    MPI_Init(&argc, &argv);
    t0 = MPI_Wtime();
    /* The default handler would abort the job on the first transport error, before anything is reported. */
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    timeout = 30.0;
    bad = argc > 4;
    if (argc > 1) {
        timeout = strtod(argv[1], &end);
        bad |= end == argv[1] || *end != '\0' || !(timeout > 0);
    }
    if (argc > 2) {
        long v = strtol(argv[2], &end, 10);

        bad |= end == argv[2] || *end != '\0' || v < 0 || v > 1000000;
        rounds = (int)v;
    }
    if (argc > 3) {
        factor = strtod(argv[3], &end);
        bad |= end == argv[3] || *end != '\0' || !(factor > 1);
    }
    if (bad) {
        if (rank == 0) {
            fprintf(stderr, "usage: %s [timeout_s > 0 [latency_rounds >= 0 [outlier_factor > 1]]]\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    memset(&mine, 0, sizeof mine);
    MPI_Get_processor_name(mine.host, &len);
    strcpy(host, mine.host);

    req = (MPI_Request *)malloc(4 * (size_t)size * sizeof *req);
    flag = (char **)malloc(4 * (size_t)size * sizeof *flag);
    idx = (int *)malloc(4 * (size_t)size * sizeof *idx);
    st = (MPI_Status *)malloc(4 * (size_t)size * sizeof *st);
    rstate = (char *)calloc(size, 1);
    sstate = (char *)calloc(size, 1);
    rbuf = (int *)malloc((size_t)size * sizeof *rbuf);
    if (!req || !flag || !idx || !st || !rstate || !sstate || !rbuf) {
        fprintf(stderr, "rank %d: out of memory\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    rstate[rank] = sstate[rank] = DONE;
    rbuf[rank] = rank;

    /* Shifted exchange keeps ~2 requests in flight on healthy links. Capping
     * each round at w bounds every rank's posting to ~timeout after its start,
     * so 2*timeout cannot flag a working link merely because its peer was
     * delayed by failures elsewhere. */
    w = size > 1 ? timeout / (size - 1) : 0.0;
    for (k = 1; k < size; ++k) {
        int src = (rank - k + size) % size, dst = (rank + k) % size;
        double tr = MPI_Wtime();

        post(0, &rbuf[src], 1, MPI_INT, src, TAG_PROBE, &rstate[src]);
        post(1, &rank, 1, MPI_INT, dst, TAG_PROBE, &sstate[dst]);
        while ((rstate[src] == PENDING || sstate[dst] == PENDING) && MPI_Wtime() - tr < w) {
            progress();
        }
    }
    while (nact > 0 && MPI_Wtime() - t0 < 2 * timeout) {
        progress();
    }

    for (j = 0; j < size; ++j) {
        if (rstate[j] == DONE && rbuf[j] != j) {
            rstate[j] = FAILED;
        }
    }

    /* Only a completed receive proves a link: eager sends may complete without delivery. */
    nfail = count_not_done(rstate, size);
    mine.nfail = nfail;
    j = 0;
    while (mine.nranges < MAXRANGES && next_run(rstate, size, &j, &lo, &hi)) {
        mine.range[mine.nranges][0] = lo;
        mine.range[mine.nranges][1] = hi;
        ++mine.nranges;
    }
    /* stdout is forwarded by the launcher, not the MPI fabric, so this hopefully survives broken links. */
    if (nfail > 0) {
        printf("rank %d (%s): no message from %d rank(s): ", rank, mine.host, nfail);
        print_runs(rstate, size);
        printf("\n");
    }
    nsend = count_not_done(sstate, size);
    if (nsend > 0) {
        printf("rank %d (%s): %d send(s) not completed: ", rank, mine.host, nsend);
        print_runs(sstate, size);
        printf("\n");
    }
    fflush(stdout);

    if (rank != 0) {
        char go = PENDING;
        int gobuf;

        post(1, &mine, (int)sizeof mine, MPI_BYTE, 0, TAG_REPORT, &rep);
        while (nact > 0 && MPI_Wtime() - t0 < 3 * timeout) {
            progress();
        }
        if (rep != DONE) {
            printf("rank %d (%s): report to rank 0 not completed\n", rank, mine.host);
            fflush(stdout);
        }
        if (!(nact == 0 && nfail == 0 && rep == DONE && count_not_done(sstate, size) == 0)) {
            /* Finalizing with requests pending is erroneous. Rank 0 normally
             * aborts the job first; this is the fallback if it can't. */
            while (MPI_Wtime() - t0 < 4 * timeout) {
                progress();
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (rounds == 0) {
            MPI_Finalize();
            return 0;
        }
        /* Only rank 0 knows whether every link works; latency rounds must not start otherwise. */
        post(0, &gobuf, 1, MPI_INT, 0, TAG_GO, &go);
        while (nact > 0 && MPI_Wtime() - t0 < 4 * timeout) {
            progress();
        }
        if (go != DONE) {
            printf("rank %d (%s): no go-ahead from rank 0 for latency rounds\n", rank, mine.host);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    } else {
        report_t *all = (report_t *)calloc(size, sizeof *all);
        char *repstate = (char *)calloc(size, 1), *unheard = (char *)calloc(size, 1), *gostate;
        int *order = (int *)malloc((size_t)size * sizeof *order);
        int nmiss, nknown = 0, nhosts = 0, ntrunc = 0, a, b, i;
        long npairs = 0;
        double elapsed, tg;

        if (!all || !repstate || !unheard || !order) {
            fprintf(stderr, "rank 0: out of memory\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        all[0] = mine;
        repstate[0] = DONE;
        for (j = 1; j < size; ++j) {
            post(0, &all[j], (int)sizeof all[j], MPI_BYTE, j, TAG_REPORT, &repstate[j]);
        }
        while (nact > 0 && MPI_Wtime() - t0 < 3 * timeout) {
            progress();
        }
        elapsed = MPI_Wtime() - t0;

        nmiss = count_not_done(repstate, size);
        for (j = 0; j < size; ++j) {
            int covered = 0;

            if (repstate[j] != DONE) {
                continue;
            }
            order[nknown++] = j;
            npairs += all[j].nfail;
            for (k = 0; k < all[j].nranges; ++k) {
                for (i = all[j].range[k][0]; i <= all[j].range[k][1]; ++i) {
                    unheard[i] = 1;
                    ++covered;
                }
            }
            if (covered < all[j].nfail) {
                ++ntrunc;
            }
        }
        g_rep = all;
        qsort(order, nknown, sizeof *order, by_host);
        for (k = 0; k < nknown; ++k) {
            if (k == 0 || strcmp(all[order[k]].host, all[order[k - 1]].host) != 0) {
                ++nhosts;
            }
        }

        if (npairs > 0 || nmiss > 0) {
            printf("connectivity FAILED (%.1f s): %ld directed rank pair(s) without delivery, %d report(s) missing\n", elapsed, npairs, nmiss);
            for (a = 0; a < nknown; a = b) {
                int n = 0, nf = 0, nu = 0;

                for (b = a; b < nknown && strcmp(all[order[b]].host, all[order[a]].host) == 0; ++b) {
                    ++n;
                    nf += all[order[b]].nfail > 0;
                    nu += unheard[order[b]];
                }
                if (nf > 0 || nu > 0) {
                    printf("  %s: %d rank(s), %d missed messages from some peer, %d not heard by some peer\n", all[order[a]].host, n, nf, nu);
                }
            }
            if (nmiss > 0) {
                printf("  no report, host unknown, from rank(s): ");
                print_runs(repstate, size);
                printf("\n");
            }
            if (ntrunc > 0) {
                printf("  %d report(s) truncated; the per-rank lines have complete lists\n", ntrunc);
            }
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        printf("connectivity OK: each of %d rank(s) on %d host(s) received from every other rank (%.2f s)\n", size, nhosts, elapsed);
        fflush(stdout);
        if (rounds == 0) {
            MPI_Finalize();
            return 0;
        }
        gostate = (char *)calloc(size, 1);
        if (!gostate) {
            fprintf(stderr, "rank 0: out of memory\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (j = 1; j < size; ++j) {
            post(1, &rank, 1, MPI_INT, j, TAG_GO, &gostate[j]);
        }
        tg = MPI_Wtime();
        while (nact > 0 && MPI_Wtime() - tg < timeout) {
            progress();
        }
    }

    nmine = size - 1 - rank;
    samp = (float *)malloc(((size_t)nmine * rounds + 1) * sizeof *samp);
    pstep = (int *)calloc((size_t)nmine + 1, sizeof *pstep);
    dead = (char *)calloc(size, 1);
    inbuf = (int *)malloc((size_t)size * sizeof *inbuf);
    mystats = (pairstat_t *)calloc((size_t)nmine + 1, sizeof *mystats);
    ev = (event_t *)calloc(MAXEV, sizeof *ev);
    if (!samp || !pstep || !dead || !inbuf || !mystats || !ev) {
        fprintf(stderr, "rank %d: out of memory for latency rounds\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    /* Negative marks a round without a sample. */
    for (z = 0; z < (size_t)nmine * rounds; ++z) {
        samp[z] = -1.0f;
    }
    tl = MPI_Wtime();
    hdr[0] = latency(rounds, samp, pstep, dead, inbuf, &hdr[1]);
    tl_end = MPI_Wtime();
    hdr[2] = pair_stats(rounds, factor, samp, pstep, dead, mystats, nmine, ev);

    if (rank != 0) {
        char hflag = PENDING, eflag = PENDING, sflag = DONE;

        post(1, hdr, 3, MPI_INT, 0, TAG_HDR, &hflag);
        post(1, ev, MAXEV * (int)sizeof *ev, MPI_BYTE, 0, TAG_EV, &eflag);
        if (nmine > 0) {
            sflag = PENDING;
            post(1, mystats, nmine * (int)sizeof *mystats, MPI_BYTE, 0, TAG_STATS, &sflag);
        }
        while (nact > 0 && MPI_Wtime() - tl_end < 3 * timeout) {
            progress();
        }
        if (nact == 0 && hflag == DONE && eflag == DONE && sflag == DONE && hdr[0] == 0 && !hdr[1]) {
            MPI_Finalize();
            return 0;
        }
        while (MPI_Wtime() - tl_end < 5 * timeout) {
            progress();
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    } else {
        long nps = (long)size * (size - 1) / 2, off = nmine, nuse = 0, nev = 0;
        pairstat_t *ps = (pairstat_t *)calloc((size_t)nps + 1, sizeof *ps), *use = (pairstat_t *)calloc((size_t)nps + 1, sizeof *use);
        event_t *evall = (event_t *)calloc((size_t)size * MAXEV, sizeof *evall), *evuse = (event_t *)calloc((size_t)size * MAXEV, sizeof *evuse);
        int *hdrs = (int *)calloc(3 * (size_t)size, sizeof *hdrs), nmiss = 0, nstall = 0, ngave = 0, ntrunc = 0, outliers;
        char *hst = (char *)calloc(size, 1), *est = (char *)calloc(size, 1), *sst = (char *)calloc(size, 1);

        if (!ps || !use || !evall || !evuse || !hdrs || !hst || !est || !sst) {
            fprintf(stderr, "rank 0: out of memory for latency statistics\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        memcpy(ps, mystats, (size_t)nmine * sizeof *ps);
        memcpy(evall, ev, MAXEV * sizeof *ev);
        hdrs[0] = hdr[0];
        hdrs[1] = hdr[1];
        hdrs[2] = hdr[2];
        hst[0] = est[0] = sst[0] = DONE;
        for (j = 1; j < size; ++j) {
            int nj = size - 1 - j;

            post(0, &hdrs[3 * j], 3, MPI_INT, j, TAG_HDR, &hst[j]);
            post(0, &evall[(size_t)j * MAXEV], MAXEV * (int)sizeof *evall, MPI_BYTE, j, TAG_EV, &est[j]);
            if (nj > 0) {
                post(0, &ps[off], nj * (int)sizeof *ps, MPI_BYTE, j, TAG_STATS, &sst[j]);
            } else {
                sst[j] = DONE;
            }
            off += nj;
        }
        while (nact > 0 && MPI_Wtime() - tl_end < 3 * timeout) {
            progress();
        }
        /* Copy out only completed receives: a pending receive's buffer still belongs to MPI, and RDMA may write it at any time. */
        for (off = 0, j = 0; j < size; off += size - 1 - j, ++j) {
            int nj = size - 1 - j;
            long nsp = 0;

            if (hst[j] != DONE || est[j] != DONE || sst[j] != DONE) {
                ++nmiss;
                continue;
            }
            nstall += hdrs[3 * j];
            ngave += hdrs[3 * j + 1] != 0;
            memcpy(&use[nuse], &ps[off], (size_t)nj * sizeof *use);
            for (k = 0; k < nj; ++k) {
                nsp += ps[off + k].nspike;
            }
            nuse += nj;
            memcpy(&evuse[nev], &evall[(size_t)j * MAXEV], (size_t)hdrs[3 * j + 2] * sizeof *evuse);
            nev += hdrs[3 * j + 2];
            ntrunc += nsp > hdrs[3 * j + 2];
        }
        outliers = analyze(use, nuse, evuse, nev, ntrunc, rounds, factor, tl_end - tl);
        if (nmiss > 0 || nstall > 0 || ngave > 0) {
            printf("latency INCOMPLETE: %d stalled exchange(s), %d rank(s) gave up, %d rank(s) sent no statistics\n", nstall, ngave, nmiss);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        MPI_Finalize();
        return outliers ? 2 : 0;
    }
    return 1;
}
