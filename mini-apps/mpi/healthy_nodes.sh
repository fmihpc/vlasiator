# shellcheck shell=bash
# Source this from a Slurm batch script (submitted with --no-kill), then call:  find_healthy_nodes MIN MAX
# It returns 0 and sets HEALTHY_NODES (comma-separated, for srun --nodelist) and HEALTHY_COUNT, or returns 1 with the reason on stderr.
#
# Optional settings, read when the function is called:
#   HEALTH_PROBE_TIMEOUT   seconds allowed per Slurm probe (default 60)
#   HEALTH_PROBE_PARALLEL  Slurm probes running at once (default 32)
#   HEALTH_CONN_PER_NODE   ranks per node in the MPI check; set it to the NICs per node to exercise each NIC (default 1)
#   HEALTH_CONN_TIMEOUT    whole seconds passed to connectivity --timeout (default 30)
#   HEALTH_CONN_ATTEMPTS   MPI checks before giving up (default 3)
#   HEALTH_CONNECTIVITY    the connectivity program (default: next to this file)
# The output of each MPI check is kept as connectivity.<jobid>.<attempt>.out in the working directory.

if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    echo "${0##*/} is meant to be sourced from a batch script: source healthy_nodes.sh; find_healthy_nodes MIN MAX" >&2
    exit 1
fi

# Resolved now, while the path used in "source" is still valid even if the caller later changes directory.
_HN_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Reads "connectivity --machine" output and prints the nodes to drop. Exit 2: no records at all; exit 3: record node names
# don't match the node list, so nothing can be concluded.
# shellcheck disable=SC2016  # the dollar signs belong to awk, not to the shell
_HN_DROP_AWK='
{
    p = index($0, "@connectivity ")
    if (p == 0) next
    n = split(substr($0, p), fld, " ")
    r = ""; node = ""; miss = ""
    for (i = 2; i <= n; ++i) {
        eq = index(fld[i], "=")
        k = substr(fld[i], 1, eq - 1); v = substr(fld[i], eq + 1)
        if (k == "rank") r = v; else if (k == "node") node = v; else if (k == "missing") miss = v
    }
    if (r == "" || node == "") next
    nodeof[r] = node; missing[r] = miss; cnt[node]++; nrec++
}
END {
    if (nrec == 0) exit 2
    nn = split(nodes, list, ",")
    known = 0
    for (i = 1; i <= nn; ++i) {
        inlist[list[i]] = 1
        if (list[i] in cnt) known++
        # Fewer records than ranks: a rank died before reporting, so this node failures are unknown.
        if (cnt[list[i]] + 0 < per) drop[list[i]] = 1
    }
    if (known == 0) exit 3
    for (r in missing) {
        if (missing[r] == "-") continue
        m = split(missing[r], part, ",")
        for (j = 1; j <= m; ++j) {
            if (index(part[j], "-")) { split(part[j], b, "-"); lo = b[1] + 0; hi = b[2] + 0 } else { lo = part[j] + 0; hi = lo }
            for (s = lo; s <= hi; ++s) {
                if (!(s in nodeof)) continue
                a = nodeof[r]; c = nodeof[s]
                # Ranks on one node failing each other: only dropping that node helps.
                if (a == c) { drop[a] = 1; continue }
                if (a > c) { t = a; a = c; c = t }
                w[a SUBSEP c]++
            }
        }
    }
    # Greedy cover of the failed node pairs: an isolated node fails with everyone, so it goes first, and
    # the nodes that merely talked to it are left without failures.
    for (;;) {
        split("", deg); split("", wt)
        for (e in w) {
            split(e, h, SUBSEP)
            if ((h[1] in drop) || (h[2] in drop)) continue
            deg[h[1]]++; deg[h[2]]++; wt[h[1]] += w[e]; wt[h[2]] += w[e]
        }
        best = ""
        for (x in deg)
            if (best == "" || deg[x] > deg[best] || (deg[x] == deg[best] && (wt[x] > wt[best] || (wt[x] == wt[best] && x > best))))
                best = x
        if (best == "") break
        drop[best] = 1
    }
    for (x in drop) if (x in inlist) print x
}'

# Prints how many distinct ranks reported from a job of size n, and a size seen that differs from n (empty if none).
# shellcheck disable=SC2016  # the dollar signs belong to awk, not to the shell
_HN_COUNT_AWK='
{
    p = index($0, "@connectivity ")
    if (p == 0) next
    m = split(substr($0, p), f, " ")
    r = ""; sz = ""
    for (i = 2; i <= m; ++i) {
        eq = index(f[i], "=")
        k = substr(f[i], 1, eq - 1); v = substr(f[i], eq + 1)
        if (k == "rank") r = v; else if (k == "size") sz = v
    }
    if (sz + 0 != n + 0) other = sz
    else if (!(r in seen)) { seen[r] = 1; good++ }
}
END { print good + 0, other }'

_hn_posint() {
    [[ $2 =~ ^[1-9][0-9]*$ ]] && return 0
    echo "find_healthy_nodes: $1 must be a positive whole number, not '$2'" >&2
    return 1
}

# shellcheck disable=SC2034  # HEALTHY_NODES and HEALTHY_COUNT are outputs for the calling script
find_healthy_nodes() {
    local min=${1:-} max=${2:-}
    local probe_timeout=${HEALTH_PROBE_TIMEOUT:-60} probe_parallel=${HEALTH_PROBE_PARALLEL:-32} per=${HEALTH_CONN_PER_NODE:-1}
    local conn_timeout=${HEALTH_CONN_TIMEOUT:-30} attempts=${HEALTH_CONN_ATTEMPTS:-3} conn=${HEALTH_CONNECTIVITY:-$_HN_DIR/connectivity}
    local attempt=0 list out rc awkrc drop_out node expected good other
    local -a all cand drop keep use

    HEALTHY_NODES=""
    HEALTHY_COUNT=0
    if [[ $# -ne 2 ]]; then
        echo "find_healthy_nodes: usage: find_healthy_nodes MIN MAX" >&2
        return 1
    fi
    _hn_posint MIN "$min" && _hn_posint MAX "$max" && _hn_posint HEALTH_PROBE_TIMEOUT "$probe_timeout" && _hn_posint HEALTH_PROBE_PARALLEL "$probe_parallel" &&
        _hn_posint HEALTH_CONN_PER_NODE "$per" && _hn_posint HEALTH_CONN_TIMEOUT "$conn_timeout" && _hn_posint HEALTH_CONN_ATTEMPTS "$attempts" || return 1
    if (( min > max )); then
        echo "find_healthy_nodes: MIN ($min) is larger than MAX ($max)" >&2
        return 1
    fi
    if [[ -z ${SLURM_JOB_NODELIST:-} ]]; then
        echo "find_healthy_nodes: SLURM_JOB_NODELIST is not set; call this inside a Slurm job" >&2
        return 1
    fi
    if [[ ! -x $conn ]]; then
        echo "find_healthy_nodes: no connectivity program at $conn; build it there or set HEALTH_CONNECTIVITY" >&2
        return 1
    fi

    # 1. Nodes on which Slurm can still start a task. A step aimed at an unresponsive node can hang until Slurm marks the node
    #    down, so every probe gets a deadline.
    mapfile -t all < <(scontrol show hostnames "$SLURM_JOB_NODELIST")
    mapfile -t cand < <(printf '%s\n' "${all[@]}" |
        xargs -P "$probe_parallel" -I{} sh -c "timeout $probe_timeout srun --nodes=1 --ntasks=1 --nodelist={} true >/dev/null 2>&1 && echo {}" |
        sort)
    echo "find_healthy_nodes: reachable through Slurm: ${#cand[@]} of ${#all[@]} node(s): ${cand[*]}"

    # 2. MPI connectivity among all of them, spares included, so a spare can replace a node the fabric can't reach.
    #    After dropping nodes the check runs again, so the nodes finally returned have passed one clean check together.
    while :; do
        attempt=$((attempt + 1))
        if (( ${#cand[@]} < min )); then
            echo "find_healthy_nodes: only ${#cand[@]} usable node(s), at least $min needed" >&2
            return 1
        fi
        list=$(IFS=,; echo "${cand[*]}")
        out=connectivity.${SLURM_JOB_ID:-local}.$attempt.out
        rc=0
        # The outer timeout only guards against a hang in MPI_Init, which the program's own deadlines can't cover.
        timeout $((5 * conn_timeout + 120)) srun --nodes="${#cand[@]}" --nodelist="$list" --ntasks-per-node="$per" "$conn" --machine --timeout="$conn_timeout" > "$out" 2>&1 || rc=$?
        grep -v '@connectivity ' "$out" || true
        if (( rc == 0 )); then
            # Separate single-rank MPI jobs also exit 0, each trivially fine, when srun and the MPI library don't agree on how
            # to start MPI. Only records from one job of the expected size show that the check tested anything.
            expected=$(( ${#cand[@]} * per ))
            read -r good other < <(awk -v n="$expected" "$_HN_COUNT_AWK" "$out")
            if [[ -n $other ]]; then
                echo "find_healthy_nodes: the MPI check ran as separate MPI jobs of size $other instead of one of size $expected, so it tested nothing." >&2
                echo "find_healthy_nodes: srun is not starting your MPI library as one job: compare 'srun --mpi=list' with how it was built, then e.g. export SLURM_MPI_TYPE=pmix" >&2
                return 1
            fi
            if (( good != expected )); then
                echo "find_healthy_nodes: the MPI check exited 0 but only $good of $expected ranks reported; see $out" >&2
                return 1
            fi
            echo "find_healthy_nodes: MPI check $attempt passed on ${#cand[@]} node(s)"
            break
        fi
        # 64 means connectivity rejected its command line, which says nothing about the nodes.
        if (( rc == 64 )); then
            echo "find_healthy_nodes: connectivity rejected its arguments; see $out" >&2
            return 1
        fi
        if (( attempt >= attempts )); then
            echo "find_healthy_nodes: MPI check still failing after $attempt attempts; see connectivity.*.out" >&2
            return 1
        fi
        awkrc=0
        drop_out=$(awk -v nodes="$list" -v per="$per" "$_HN_DROP_AWK" "$out") || awkrc=$?
        case $awkrc in
            0) ;;
            2) echo "find_healthy_nodes: MPI check $attempt left no machine-readable records in $out, so the bad node(s) can't be identified" >&2; return 1 ;;
            3) echo "find_healthy_nodes: node names in $out don't match Slurm's node list, so the bad node(s) can't be identified" >&2; return 1 ;;
            *) echo "find_healthy_nodes: could not evaluate $out" >&2; return 1 ;;
        esac
        if [[ -z $drop_out ]]; then
            echo "find_healthy_nodes: MPI check $attempt failed but no node explains it; checking again"
            continue
        fi
        mapfile -t drop < <(sort <<< "$drop_out")
        echo "find_healthy_nodes: MPI check $attempt failed; dropping: ${drop[*]}"
        keep=()
        for node in "${cand[@]}"; do
            [[ " ${drop[*]} " == *" $node "* ]] || keep+=("$node")
        done
        cand=("${keep[@]}")
    done

    use=("${cand[@]:0:$max}")
    HEALTHY_COUNT=${#use[@]}
    HEALTHY_NODES=$(IFS=,; echo "${use[*]}")
    echo "find_healthy_nodes: using $HEALTHY_COUNT node(s): ${use[*]}"
    return 0
}
