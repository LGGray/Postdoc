#!/usr/bin/env python3
"""Which genes' allelic ratio differs between two sections, calibrated against
an autosomal null.

THE QUESTION. The tile x gene table carries every gene, not just the escape
panel, so "which X-linked genes escape more in one section than the other" is a
groupby away. Pooling each gene's molecules per section and testing the CAST
fraction is the obvious thing to do, and on its own it is WRONG - the count of
"significant" genes it returns is mostly technical.

WHY THE OBVIOUS TEST FAILS, AND HOW THIS ONE IS CALIBRATED
----------------------------------------------------------
A two-proportion test assumes the molecules are independent draws. They are not:
molecules cluster within cells and tiles, ambient RNA is shared, and the two
libraries differ in depth and duplication. Run the identical test on AUTOSOMAL
genes - which are under no allelic regulation and cannot genuinely shift between
two sections - and it calls ~7% of them significant at BH q<0.05. That number is
the false-positive rate, measured rather than assumed.

So the autosomes are used as the null. Their z should be N(0,1) if the binomial
model held; the inflation factor sqrt(median(z^2)/0.4549) says how far off it
is, and it GROWS WITH DEPTH (~1.3 at 200-500 molecules, ~2.6 above 5000),
because overdispersion scales with n while binomial variance does not. The
factor is therefore measured per depth stratum and divided out of the target
chromosome's z within the matching stratum. After correction the autosomes come
back at ~0.4%, which is printed as a calibration check - if that number is not
small, the correction has not worked and nothing downstream should be believed.

This is genomic control, borrowed from GWAS. It is conservative for real signal
and it is the difference between 9 hits and 8 hits here, but more importantly it
is the difference between quoting q=1e-10 and knowing what q means.

WHAT IT DOES NOT FIX. Cell-type composition. If the two sections sample
different proportions of myocyte/fibroblast/endothelium and escape differs by
cell type, a gene's pooled ratio moves for a reason that is neither age nor
chance, and no amount of calibration sees it. Genes moving in BOTH directions
is the signature to watch for: uniform loss of silencing predicts one direction,
so a mixed set argues for composition. Report the directions.

AND THE DESIGN LIMIT. One section per age here (n = 1 animal), so every
difference this finds describes two sections, not two ages. The output carries
the sample labels, never "adult"/"aged" as a group.

LEVEL. Molecules (a1_umi/a2_umi) by default. --level reads is duplicate-
inclusive AND the duplication is allele-asymmetric, so it biases the ratio
toward the more-amplified allele; it is available for comparison and warns.

OUTPUTS  (--out-dir)
--------------------
  chrx_allelic_shift_<level>.tsv         one row per gene tested, both
                                         chromosomes' worth: raw and corrected
                                         statistics, plus the untestable genes
                                         flagged rather than dropped
  chrx_shift_tiles_<level>.tsv.gz        per-tile rows for --chrom, the input to
                                         chrx_escape_clustering.py
  chrx_allelic_shift.provenance.tsv

  sbatch ~/Postdoc/slurm/spatial_chrx_shift.slurm
"""
import argparse
import collections
import datetime
import gzip
import math
import os
import sys

MEDIAN_CHISQ_1DF = 0.4549364231195736   # median of chi2(1); the genomic-control constant


def msg(fmt, *a):
    sys.stderr.write((fmt % a if a else fmt) + "\n")
    sys.stderr.flush()


def norm_sf2(z):
    """two-sided tail of the standard normal"""
    return math.erfc(abs(z) / math.sqrt(2.0))


def bh(ps):
    """Benjamini-Hochberg q-values, in the input order"""
    m = len(ps)
    if not m:
        return []
    order = sorted(range(m), key=lambda i: ps[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        prev = min(prev, ps[i] * m / (m - rank + 1))
        q[i] = prev
    return q


def median(v):
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def open_counts(root, sample):
    """tile_gene_counts.tsv[.gz]; .gz is what the counter writes, the plain name
    is what an interrupted gunzip leaves behind - accept either, as the R scripts do."""
    f = os.path.join(root, sample, "tile_gene_counts.tsv")
    if os.path.exists(f + ".gz"):
        return gzip.open(f + ".gz", "rt"), f + ".gz"
    if os.path.exists(f):
        return open(f), f
    raise SystemExit("No tile_gene_counts.tsv[.gz] under %s" % os.path.join(root, sample))


def read_sample(root, sample, level, chrom, want_tiles):
    """One pass. Returns per-gene (a1, a2, n_tiles) and, for `chrom`, per-tile rows."""
    a1c = "a1_" + level
    a2c = "a2_" + level
    fh, path = open_counts(root, sample)
    per_gene = collections.defaultdict(lambda: [0, 0, 0])
    gene_chrom = {}
    tiles = []
    with fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(head)}
        for col in ("gene", "chrom", "trow", "tcol", a1c, a2c):
            if col not in ix:
                raise SystemExit("%s has no column %r" % (path, col))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            g = f[ix["gene"]]
            c = f[ix["chrom"]]
            a1 = int(f[ix[a1c]])
            a2 = int(f[ix[a2c]])
            rec = per_gene[g]
            rec[0] += a1
            rec[1] += a2
            rec[2] += 1
            gene_chrom[g] = c
            if want_tiles and c == chrom and (a1 + a2) > 0:
                tiles.append((g, int(f[ix["trow"]]), int(f[ix["tcol"]]), a1, a2))
    return per_gene, gene_chrom, tiles, path


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-root", required=True,
                   help="tree holding <sample>/tile_gene_counts.tsv[.gz]")
    p.add_argument("--samples", default="9w,78w",
                   help="exactly two, in the order the comparison is reported "
                        "(second minus first)")
    p.add_argument("--out-dir", default=None, help="default: --in-root")
    p.add_argument("--extra-root", default=None,
                   help="second counts tree merged over the first BY GENE, same "
                        "semantics as the R figures: a gene present here "
                        "REPLACES its rows from --in-root rather than adding")
    p.add_argument("--chrom", default="chrX", help="the chromosome under test")
    p.add_argument("--level", choices=("umi", "reads"), default="umi",
                   help="umi = molecules (default). reads is duplicate-inclusive "
                        "and allele-asymmetric; it biases the ratio")
    p.add_argument("--min-molecules", type=int, default=200,
                   help="required in BOTH samples for a gene to be tested")
    p.add_argument("--strata", default="200,500,1500,5000",
                   help="lower edges of the depth strata the inflation factor is "
                        "measured in, on min(n1, n2)")
    p.add_argument("--min-null-genes", type=int, default=30,
                   help="a stratum with fewer autosomal genes than this borrows "
                        "the nearest populated stratum's factor")
    args = p.parse_args()

    started = datetime.datetime.now()
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    if len(samples) != 2:
        raise SystemExit("--samples needs exactly two, got %r" % args.samples)
    out_dir = args.out_dir or args.in_root
    os.makedirs(out_dir, exist_ok=True)
    if args.level == "reads":
        msg("WARNING: --level reads is duplicate-inclusive and the duplication is\n"
            "         allele-asymmetric, so the ratio is biased toward the more\n"
            "         amplified allele. Molecules (umi) is the defensible level.")

    # ------------------------------------------------------------------ load
    per_gene = {}
    gene_chrom = {}
    tiles = {}
    inputs = []
    for s in samples:
        msg("Reading %s / %s", args.in_root, s)
        pg, gc, tl, path = read_sample(args.in_root, s, args.level, args.chrom, True)
        per_gene[s] = pg
        gene_chrom.update(gc)
        tiles[s] = tl
        inputs.append(path)

    if args.extra_root:
        msg("Merging extra counts from %s (by gene, replacing)", args.extra_root)
        for s in samples:
            pg, gc, tl, path = read_sample(args.extra_root, s, args.level, args.chrom, True)
            inputs.append(path)
            for g, rec in pg.items():
                per_gene[s][g] = rec          # replace, never add
                gene_chrom[g] = gc[g]
            hits = set(pg)
            if hits:
                tiles[s] = [t for t in tiles[s] if t[0] not in hits] + tl
        msg("  replaced %d gene(s): %s", len(hits), ", ".join(sorted(hits)))

    # ------------------------------------------------------------------ test
    s1, s2 = samples
    rows = []
    for g in set(per_gene[s1]) & set(per_gene[s2]):
        a1_1, a2_1, nt1 = per_gene[s1][g]
        a1_2, a2_2, nt2 = per_gene[s2][g]
        n1, n2 = a1_1 + a2_1, a1_2 + a2_2
        if n1 < args.min_molecules or n2 < args.min_molecules:
            continue
        e1, e2 = a2_1 / n1, a2_2 / n2
        pooled = (a2_1 + a2_2) / (n1 + n2)
        se = math.sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))
        r = dict(gene=g, chrom=gene_chrom[g], n_1=n1, n_2=n2, nt_1=nt1, nt_2=nt2,
                 a2_1=a2_1, a2_2=a2_2, esc_1=e1, esc_2=e2, d_esc=e2 - e1,
                 depth=min(n1, n2))
        if se == 0:
            # No minor-allele molecule anywhere: fully silenced (or fully CAST).
            # Kept in the table and flagged, because "no signal" and "no shift"
            # are different statements and dropping them hides the former.
            r.update(z_raw=float("nan"), p_raw=float("nan"), testable=0)
        else:
            z = (e2 - e1) / se
            r.update(z_raw=z, p_raw=norm_sf2(z), testable=1)
        rows.append(r)

    tested = [r for r in rows if r["testable"]]
    auto = [r for r in tested if r["chrom"] not in (args.chrom, "chrY", "chrM")]
    targ = [r for r in tested if r["chrom"] == args.chrom]
    if not auto:
        raise SystemExit("No autosomal genes passed --min-molecules; cannot calibrate")

    # ------------------------------------------------- inflation, per stratum
    edges = [int(x) for x in args.strata.split(",")] + [10 ** 18]
    strata = list(zip(edges[:-1], edges[1:]))
    lam = {}
    msg("")
    msg("Inflation of the test, measured on autosomes (1.0 = binomial holds):")
    for lo, hi in strata:
        zs = [r["z_raw"] for r in auto if lo <= r["depth"] < hi]
        if len(zs) < args.min_null_genes:
            msg("  depth %-12s %5d autosomal genes   too few, borrows a neighbour",
                "%d-%d" % (lo, hi) if hi < 10 ** 18 else "%d+" % lo, len(zs))
            continue
        lam[(lo, hi)] = math.sqrt(median([z * z for z in zs]) / MEDIAN_CHISQ_1DF)
        msg("  depth %-12s %5d autosomal genes   inflation %.2f",
            "%d-%d" % (lo, hi) if hi < 10 ** 18 else "%d+" % lo, len(zs), lam[(lo, hi)])
    if not lam:
        raise SystemExit("No depth stratum had >= %d autosomal genes" % args.min_null_genes)

    def inflation(depth):
        for (lo, hi), v in lam.items():
            if lo <= depth < hi:
                return v
        # nearest populated stratum by lower edge
        return lam[min(lam, key=lambda k: abs(k[0] - depth))]

    for r in tested:
        r["inflation"] = inflation(r["depth"])
        r["z_corr"] = r["z_raw"] / r["inflation"]
        r["p_corr"] = norm_sf2(r["z_corr"])

    # BH within each set separately: the target chromosome is the experiment,
    # the autosomes are the control, and pooling them would let the control's
    # size set the experiment's threshold.
    for sel in (targ, auto):
        for key in ("p_raw", "p_corr"):
            qs = bh([r[key] for r in sel])
            for r, q in zip(sel, qs):
                r[key.replace("p_", "q_")] = q

    n_auto_sig = sum(1 for r in auto if r["q_corr"] < 0.05)
    n_auto_raw = sum(1 for r in auto if r["q_raw"] < 0.05)
    sig = sorted([r for r in targ if r["q_corr"] < 0.05], key=lambda r: r["p_corr"])

    msg("")
    msg("%s: %d genes at >= %d molecules in both, %d testable (%d have no minor "
        "allele at all)", args.chrom, len(targ) + sum(1 for r in rows
        if r["chrom"] == args.chrom and not r["testable"]), args.min_molecules,
        len(targ), sum(1 for r in rows if r["chrom"] == args.chrom and not r["testable"]))
    msg("%-14s %9s %9s %9s %8s %8s %9s", "gene", "esc " + s1, "esc " + s2,
        "d_esc", "z_raw", "z_corr", "q_corr")
    for r in sig:
        msg("%-14s %9.3f %9.3f %+9.3f %8.1f %8.1f %9.1e", r["gene"], r["esc_1"],
            r["esc_2"], r["d_esc"], r["z_raw"], r["z_corr"], r["q_corr"])
    up = sum(1 for r in sig if r["d_esc"] > 0)
    msg("  %d significant: %d higher in %s, %d lower", len(sig), up, s2, len(sig) - up)
    if len(sig) and up not in (0, len(sig)):
        msg("  BOTH DIRECTIONS - read the docstring: uniform loss of silencing")
        msg("  predicts one direction, so this pattern is what cell-type")
        msg("  composition between two sections looks like.")
    msg("  calibration check: autosomes %d/%d (%.1f%%) before correction -> "
        "%d/%d (%.1f%%) after", n_auto_raw, len(auto), 100 * n_auto_raw / len(auto),
        n_auto_sig, len(auto), 100 * n_auto_sig / len(auto))
    if len(auto) and 100 * n_auto_sig / len(auto) > 2.0:
        msg("  WARNING: the corrected autosomal rate is still high. The")
        msg("  correction has not absorbed the technical difference between")
        msg("  these two libraries; treat every hit below as unproven.")

    # --------------------------------------------------------------- outputs
    f_tsv = os.path.join(out_dir, "chrx_allelic_shift_%s.tsv" % args.level)
    cols = ["gene", "chrom", "testable", "n_%s" % s1, "n_%s" % s2,
            "tiles_%s" % s1, "tiles_%s" % s2, "minor_%s" % s1, "minor_%s" % s2,
            "esc_%s" % s1, "esc_%s" % s2, "d_esc", "depth",
            "z_raw", "p_raw", "q_raw", "inflation", "z_corr", "p_corr", "q_corr"]
    key_of = {"n_%s" % s1: "n_1", "n_%s" % s2: "n_2", "tiles_%s" % s1: "nt_1",
              "tiles_%s" % s2: "nt_2", "minor_%s" % s1: "a2_1",
              "minor_%s" % s2: "a2_2", "esc_%s" % s1: "esc_1",
              "esc_%s" % s2: "esc_2"}
    with open(f_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["chrom"] != args.chrom,
                                             r.get("p_corr", 1.0))):
            out = []
            for c in cols:
                v = r.get(key_of.get(c, c), "")
                if isinstance(v, float):
                    out.append("NA" if v != v else "%.6g" % v)
                else:
                    out.append(str(v))
            fh.write("\t".join(out) + "\n")
    msg("Wrote %s (%d rows)", f_tsv, len(rows))

    f_tiles = os.path.join(out_dir, "chrx_shift_tiles_%s.tsv.gz" % args.level)
    n_tile_rows = 0
    with gzip.open(f_tiles, "wt") as fh:
        fh.write("sample\tgene\ttrow\ttcol\ta1\ta2\tn\n")
        for s in samples:
            for g, r_, c_, a1, a2 in tiles[s]:
                fh.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % (s, g, r_, c_, a1, a2, a1 + a2))
                n_tile_rows += 1
    msg("Wrote %s (%d rows)", f_tiles, n_tile_rows)

    prov = os.path.join(out_dir, "chrx_allelic_shift.provenance.tsv")
    with open(prov, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in (
            ("argv", " ".join(sys.argv)),
            ("run_at", started.strftime("%Y-%m-%d %H:%M:%S")),
            ("elapsed_s", "%.0f" % (datetime.datetime.now() - started).total_seconds()),
            ("in_root", args.in_root),
            ("extra_root", args.extra_root or "(none)"),
            ("inputs", ";".join(inputs)),
            ("samples", ",".join(samples)),
            ("chrom", args.chrom),
            ("level", args.level),
            ("min_molecules", str(args.min_molecules)),
            ("strata", args.strata),
            ("inflation_by_stratum",
             ";".join("%d-%s:%.3f" % (lo, "inf" if hi > 10 ** 17 else str(hi), v)
                      for (lo, hi), v in sorted(lam.items()))),
            ("genes_target_tested", str(len(targ))),
            ("genes_target_significant", str(len(sig))),
            ("genes_target_up", str(up)),
            ("genes_target_down", str(len(sig) - up)),
            ("genes_autosomal_null", str(len(auto))),
            ("autosomal_sig_frac_raw", "%.4f" % (n_auto_raw / len(auto))),
            ("autosomal_sig_frac_corrected", "%.4f" % (n_auto_sig / len(auto))),
            ("significant_genes", ",".join(r["gene"] for r in sig)),
        ):
            fh.write("%s\t%s\n" % (k, v))
    msg("Wrote %s", prov)


if __name__ == "__main__":
    main()
