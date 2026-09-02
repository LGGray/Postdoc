#!/usr/bin/env python3
"""Is a gene's escape spatially patchy or scattered, per section?

THE QUESTION, AND THE TWO WAYS TO GET IT WRONG
----------------------------------------------
"Escape is clustered" means: of the tiles where this gene is seen at all, the
ones carrying an inactive-X molecule sit next to each other more often than
chance. Two confounds make the naive version of that test answer a different
question.

1. THE EXPRESSION FOOTPRINT. Tissue is structured, so every gene's detected
   tiles are already clustered - Smpx is cardiomyocyte-restricted, Bgn tracks
   ECM. Comparing escape-positive tiles against the whole slide measures where
   the GENE is, not where the ESCAPE is. So the permutation null shuffles escape
   labels only among the tiles where that gene was detected in that section.

2. TILE DEPTH. A tile holding ten molecules is far likelier to catch one CAST
   molecule than a tile holding one - for Smpx the escape-positive rate runs
   from 7% at one molecule to 57% at ten or more - and depth is itself spatially
   clustered by tissue density and section thickness. So a depth-blind null
   reports clustering that is really depth clustering. --stratify-depth (on by
   default) shuffles WITHIN depth strata, which holds the depth-escape
   relationship fixed and asks whether anything spatial is left over. It matters:
   it took Bgn at 78w from z = +1.9 to +0.1, and it halves Smpx at 78w.

Both nulls are reported, unstratified first, so the difference between them is
visible instead of being a choice made off-screen.

STATISTIC. A join count: adjacent (8-neighbour) tile pairs where BOTH are
escape-positive, against a permutation null. Reported as z and an empirical p.
Join counts rather than Moran's I because per-tile escape is overwhelmingly
one molecule - a continuous ratio per tile is 0 or 1 and its autocorrelation is
not interpretable - so presence/absence is the honest resolution.

READ THE POWER COLUMN BEFORE READING THE p. A gene with 43 escape-positive
tiles among 518 detected has an EXPECTED join count near 1.8; observing 0 gives
p = 1.0 and means nothing. Those rows say "underpowered", and underpowered is
not the same as dispersed - it is no evidence either way. Only genes whose
expected join count is comfortably above zero can support a claim in either
direction.

WHAT CLUSTERING DOES NOT DISTINGUISH. Clonal XCI patches from cell-type
domains. If escape is cell-type specific and the cell type is regional, escape
clusters with no clonal patch anywhere. Separating them needs the cell-type
annotation (spatial_008um_seurat_object_banksy.RDS) and a within-domain test.

Input is chrx_shift_tiles_<level>.tsv.gz from chrx_allelic_shift.py.

OUTPUTS  (--out-dir)
--------------------
  chrx_escape_clustering.tsv               one row per gene x sample x test
  chrx_escape_clustering.provenance.tsv
"""
import argparse
import collections
import csv
import datetime
import gzip
import math
import os
import random
import sys

# Half the 8-neighbourhood, so each adjacent pair is visited exactly once.
NEIGHBOURS = ((0, 1), (1, -1), (1, 0), (1, 1))


def msg(fmt, *a):
    sys.stderr.write((fmt % a if a else fmt) + "\n")
    sys.stderr.flush()


def edge_list(positions):
    """indices of adjacent pairs among `positions` (list of (row, col))"""
    ix = {p: i for i, p in enumerate(positions)}
    out = []
    for (r, c) in positions:
        i = ix[(r, c)]
        for dr, dc in NEIGHBOURS:
            j = ix.get((r + dr, c + dc))
            if j is not None:
                out.append((i, j))
    return out


def join_count(edges, labels):
    return sum(1 for i, j in edges if labels[i] and labels[j])


def permute(positions, labels, strata, nperm, rng):
    """Observed join count against a null that shuffles labels among positions,
    within strata when given. Returns (obs, mean, sd, z, p, n_edges)."""
    edges = edge_list(positions)
    x = list(labels)
    obs = join_count(edges, x)
    if strata is None:
        groups = {0: list(range(len(x)))}
    else:
        groups = collections.defaultdict(list)
        for i, s in enumerate(strata):
            groups[s].append(i)
    held = {k: [x[i] for i in idx] for k, idx in groups.items()}
    null = []
    y = list(x)
    for _ in range(nperm):
        for k, idx in groups.items():
            vals = held[k][:]
            rng.shuffle(vals)
            for i, v in zip(idx, vals):
                y[i] = v
        null.append(join_count(edges, y))
    mu = sum(null) / len(null)
    var = sum((v - mu) ** 2 for v in null) / len(null)
    sd = math.sqrt(var)
    z = (obs - mu) / sd if sd > 0 else 0.0
    p = (sum(1 for v in null if v >= obs) + 1) / (nperm + 1)
    return obs, mu, sd, z, p, len(edges)


def depth_bin(n, edges):
    """exact count while it is small, then coarsening bins - the escape rate
    changes fastest between one and a few molecules, which is where the
    stratification has to be finest"""
    for e in edges:
        if n < e:
            return n if n < edges[0] else e
    return edges[-1]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tiles", required=True,
                   help="chrx_shift_tiles_<level>.tsv.gz from chrx_allelic_shift.py")
    p.add_argument("--shift-table", default=None,
                   help="chrx_allelic_shift_<level>.tsv; genes are taken from it "
                        "at --q-cutoff unless --genes is given")
    p.add_argument("--genes", default=None, help="comma-separated, overrides --shift-table")
    p.add_argument("--q-cutoff", type=float, default=0.05,
                   help="q_corr below which a gene is carried forward (default 0.05)")
    p.add_argument("--out-dir", default=None, help="default: alongside --tiles")
    p.add_argument("--nperm", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--min-labelled", type=int, default=10,
                   help="a test needs at least this many tiles on BOTH sides of "
                        "the label, else it is reported as underpowered")
    p.add_argument("--depth-strata", default="5,10,20",
                   help="coarsening edges; below the first edge the exact "
                        "molecule count is its own stratum")
    p.add_argument("--no-footprint", action="store_true",
                   help="skip the detection-footprint test")
    args = p.parse_args()

    started = datetime.datetime.now()
    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.tiles))
    os.makedirs(out_dir, exist_ok=True)
    rng = random.Random(args.seed)
    strata_edges = [int(x) for x in args.depth_strata.split(",")]

    # ------------------------------------------------------------------ genes
    if args.genes:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
        picked = "--genes"
    else:
        if not args.shift_table:
            raise SystemExit("give --genes or --shift-table")
        genes = []
        with open(args.shift_table) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("q_corr", "NA") in ("NA", ""):
                    continue
                if float(row["q_corr"]) < args.q_cutoff:
                    genes.append(row["gene"])
        picked = "%s at q_corr < %g" % (os.path.basename(args.shift_table), args.q_cutoff)
    if not genes:
        raise SystemExit("no genes selected (%s)" % picked)
    msg("Testing %d gene(s) from %s: %s", len(genes), picked, ", ".join(genes))

    # ------------------------------------------------------------------ tiles
    want = set(genes)
    per = collections.defaultdict(list)     # (sample, gene) -> [(r, c, a2, n)]
    universe = collections.defaultdict(set)  # sample -> every covered tile
    op = gzip.open if args.tiles.endswith(".gz") else open
    with op(args.tiles, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s, g = row["sample"], row["gene"]
            r_, c_ = int(row["trow"]), int(row["tcol"])
            n, a2 = int(row["n"]), int(row["a2"])
            universe[s].add((r_, c_))
            if g in want and n > 0:
                per[(s, g)].append((r_, c_, a2, n))
    samples = sorted(universe, key=lambda s: (s != "9w", s))
    msg("Loaded %s: %s", args.tiles,
        ", ".join("%s %d covered tiles" % (s, len(universe[s])) for s in samples))

    out = []

    def run(gene, sample, positions, labels, strata, test, note=""):
        k = sum(labels)
        row = dict(gene=gene, sample=sample, test=test, n_positions=len(positions),
                   n_labelled=k, note=note)
        if k < args.min_labelled or len(positions) - k < args.min_labelled:
            row.update(obs_jc="NA", exp_jc="NA", sd_jc="NA", z="NA", p="NA",
                       n_edges="NA", underpowered=1)
            out.append(row)
            return row
        obs, mu, sd, z, pv, ne = permute(positions, labels, strata, args.nperm, rng)
        row.update(obs_jc=obs, exp_jc="%.2f" % mu, sd_jc="%.2f" % sd,
                   z="%.3f" % z, p="%.5f" % pv, n_edges=ne, underpowered=0)
        # An expected count this small cannot support a claim either way; the
        # p-value is real but uninformative, so say so in the row itself.
        if mu < 2.0:
            row["note"] = (note + ";" if note else "") + "exp_jc<2, uninformative"
        out.append(row)
        return row

    hdr = "%-14s %-5s %-24s %7s %7s %7s %9s %8s %8s"
    msg("")
    msg(hdr, "gene", "age", "test", "tiles", "labelled", "obsJC", "expJC", "z", "p")
    for g in genes:
        for s in samples:
            v = per.get((s, g), [])
            if not v:
                continue
            pos = [(r_, c_) for r_, c_, a2, n in v]
            esc = [1 if a2 > 0 else 0 for r_, c_, a2, n in v]
            dep = [depth_bin(n, strata_edges) for r_, c_, a2, n in v]
            for test, strata in (("escape_unstratified", None),
                                 ("escape_depth_stratified", dep)):
                r = run(g, s, pos, esc, strata, test)
                msg(hdr, g, s, test, r["n_positions"], r["n_labelled"],
                    r["obs_jc"], r["exp_jc"], r["z"], r["p"])
            if not args.no_footprint:
                # Context, not the answer: does the gene's own expression cluster?
                uni = sorted(universe[s])
                det = {(r_, c_) for r_, c_, a2, n in v}
                lab = [1 if q in det else 0 for q in uni]
                r = run(g, s, uni, lab, None, "detection_footprint")
                msg(hdr, g, s, "detection_footprint", r["n_positions"],
                    r["n_labelled"], r["obs_jc"], r["exp_jc"], r["z"], r["p"])

    cols = ["gene", "sample", "test", "n_positions", "n_labelled", "n_edges",
            "obs_jc", "exp_jc", "sd_jc", "z", "p", "underpowered", "note"]
    f_tsv = os.path.join(out_dir, "chrx_escape_clustering.tsv")
    with open(f_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    msg("")
    msg("Wrote %s (%d rows)", f_tsv, len(out))

    prov = os.path.join(out_dir, "chrx_escape_clustering.provenance.tsv")
    with open(prov, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in (
            ("argv", " ".join(sys.argv)),
            ("run_at", started.strftime("%Y-%m-%d %H:%M:%S")),
            ("elapsed_s", "%.0f" % (datetime.datetime.now() - started).total_seconds()),
            ("tiles", args.tiles),
            ("shift_table", args.shift_table or "(none)"),
            ("genes", ",".join(genes)),
            ("gene_source", picked),
            ("nperm", str(args.nperm)),
            ("seed", str(args.seed)),
            ("depth_strata", args.depth_strata),
            ("min_labelled", str(args.min_labelled)),
            ("rows", str(len(out))),
            ("underpowered_rows", str(sum(1 for r in out if r["underpowered"]))),
        ):
            fh.write("%s\t%s\n" % (k, v))
    msg("Wrote %s", prov)


if __name__ == "__main__":
    main()
