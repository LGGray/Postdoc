#!/usr/bin/env python3
"""Per-2um-bin allelic counts from a Visium HD BAM - one pass, no BAM splitting.

This is step 0/1 of choosing a tile size for the spatial allele-specific
analysis. The point is that you do NOT need to run sinto + Allelome.PRO2 once
per candidate tile size. Allelic counts at the 2um grid are the finest thing
any tiling can be built from, so you count once here and then aggregate to any
nested tile size in silico (ase_tile_sweep.R). sinto only gets run later, once,
at the size the sweep picks.

Output is one row per 2um bin that carries at least one informative UMI:

    barcode  array_row  array_col  x_ref  x_alt  a_ref  a_alt

x_* are chrX, a_* are the autosomal control set. ref = C57BL/6NJ (the allele in
column REF of the SNP bed), alt = CAST/EiJ.

Why an autosomal control. We are mapping F1 B6xCAST reads to the standard 10x
B6 reference, so CAST-allele reads are slightly less mappable and every ratio
sits a bit below its true value. That bias is common to all bins, so it barely
touches the tile-to-tile contrast the sweep measures - but running the identical
pipeline on autosomes gives an empirical null that already contains the bias,
the mapping artefacts and any technical spatial structure. chrX structure is
only interesting where it exceeds that null.

Counting unit is the UMI, not the read. PCR duplicates are perfectly correlated
in allele, so counting reads inflates n and makes every downstream binomial
anticonservative. Reads sharing (CB, UB) are collapsed and the allele is decided
by majority vote across whatever SNPs that UMI's reads happen to cover; ties are
discarded.

Run: see slurm/spatial_ase_sweep.slurm
"""

import argparse
import bisect
import collections
import csv
import gzip
import os
import sys

import pysam


# ---------------------------------------------------------------- SNP loading

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def sniff_snp_layout(path, explicit=None):
    """Work out which columns hold the ref and alt bases.

    The lab's SNP beds have been written by a few different scripts over the
    years (see OCM_heart/core_escape_SNPs.R) and the allele columns are not in a
    fixed place. Rather than hardcode a guess that fails silently, detect the
    layout and report it, and let --snp-cols override.
    """
    if explicit:
        ref_i, alt_i = (int(x) - 1 for x in explicit.split(","))
        return ref_i, alt_i, None

    bases = set("ACGT")
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            break
        else:
            raise SystemExit(f"No data lines in {path}")

    # Two separate single-base columns, e.g. chr pos pos+1 A G
    singles = [i for i, v in enumerate(f) if i >= 3 and v.upper() in bases]
    if len(singles) >= 2:
        return singles[0], singles[1], None

    # One packed column with a separator, e.g. chr pos pos+1 A/G (also A|G, A>G)
    for i, v in enumerate(f):
        if i < 3:
            continue
        for sep in ("/", "|", ">", ":"):
            if sep in v:
                a, b = v.split(sep)[:2]
                if a.upper() in bases and b.upper() in bases:
                    return i, i, sep

    # Allelome.PRO's own BED6: the name field is the two bases run together,
    # e.g. chr1 3057935 3057936 AG 0 . - first character REF, second ALT.
    for i, v in enumerate(f):
        if i < 3:
            continue
        u = v.upper()
        if len(u) == 2 and u[0] in bases and u[1] in bases:
            return i, i, ""

    raise SystemExit(
        f"Cannot find ref/alt base columns in {path}\nFirst line: {f}\n"
        "Pass --snp-cols REF,ALT (1-based) to set them explicitly."
    )


def load_snps(path, chroms, pos_col, ref_i, alt_i, packed_sep, offset, swap):
    """chrom -> (sorted position array, {pos: (ref, alt)}), 0-based positions."""
    want = set(chroms)
    tables = {c: {} for c in chroms}
    skipped = 0
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[0]
            if chrom not in want:
                continue
            try:
                pos = int(f[pos_col]) + offset
            except (ValueError, IndexError):
                skipped += 1
                continue
            if packed_sep == "":
                packed = f[ref_i].upper()
                ref, alt = packed[0], packed[1]
            elif packed_sep is not None:
                parts = f[ref_i].split(packed_sep)
                ref, alt = parts[0].upper(), parts[1].upper()
            else:
                ref, alt = f[ref_i].upper(), f[alt_i].upper()
            if swap:
                ref, alt = alt, ref
            if len(ref) != 1 or len(alt) != 1 or ref == alt:
                skipped += 1
                continue
            tables[chrom][pos] = (ref, alt)
    out = {}
    for c, d in tables.items():
        out[c] = (sorted(d), d)
    return out, skipped


# ------------------------------------------------------------ bin coordinates

def load_positions(spatial_dir):
    """barcode -> (array_row, array_col), tissue-covered bins only."""
    pq = os.path.join(spatial_dir, "tissue_positions.parquet")
    if os.path.exists(pq):
        try:
            import pyarrow.parquet as pyparquet
        except ImportError:
            raise SystemExit(
                f"Found {pq} but pyarrow is not installed in this env.\n"
                "Either `conda install -c conda-forge pyarrow`, or convert once with:\n"
                "  Rscript -e 'arrow::write_tsv_arrow(arrow::read_parquet(\"%s\"), "
                "\"%s\")'" % (pq, os.path.join(spatial_dir, "tissue_positions.tsv"))
            )
        t = pyparquet.read_table(pq).to_pydict()
        cols = {k.lower(): k for k in t}
        bc = t[cols["barcode"]]
        it = t[cols["in_tissue"]]
        rr = t[cols["array_row"]]
        cc = t[cols["array_col"]]
        return {b: (r, c) for b, i, r, c in zip(bc, it, rr, cc) if int(i) == 1}

    for name in ("tissue_positions.tsv", "tissue_positions.csv",
                 "tissue_positions_list.csv"):
        p = os.path.join(spatial_dir, name)
        if not os.path.exists(p):
            continue
        delim = "\t" if name.endswith(".tsv") else ","
        pos = {}
        with open(p) as fh:
            rdr = csv.reader(fh, delimiter=delim)
            first = next(rdr)
            has_header = not first[1].strip().lstrip("-").isdigit()
            rows = rdr if has_header else [first] + list(rdr)
            for f in rows:
                if int(f[1]) != 1:
                    continue
                pos[f[0]] = (int(f[2]), int(f[3]))
        return pos

    raise SystemExit(f"No tissue_positions.* found in {spatial_dir}")


# ------------------------------------------------------------------ main pass

def count_chrom(bam, chrom, snp_sorted, snp_map, args, stats):
    """UMI-collapsed (barcode -> [ref, alt]) for one chromosome."""
    votes = collections.defaultdict(lambda: [0, 0])   # (cb, ub) -> votes
    umi_bin = {}                                      # (cb, ub) -> cb
    n_obs = n_hit = 0

    for read in bam.fetch(chrom):
        stats["reads_seen"] += 1
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if args.drop_duplicates and read.is_duplicate:
            continue
        if read.mapping_quality < args.min_mapq:
            continue
        if not read.has_tag("CB") or not read.has_tag("UB"):
            continue
        if args.require_gene and not read.has_tag(args.gene_tag):
            continue
        cb = read.get_tag("CB")
        stats["reads_kept"] += 1

        # Most reads in a CAST cross do overlap a SNP, but skipping the ones
        # that cannot still avoids building an aligned-pairs list per read.
        lo = bisect.bisect_left(snp_sorted, read.reference_start)
        if lo >= len(snp_sorted) or snp_sorted[lo] >= read.reference_end:
            continue

        key = (cb, read.get_tag("UB"))
        seq = read.query_sequence
        qual = read.query_qualities
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            allele = snp_map.get(rpos)
            if allele is None:
                continue
            if qual is not None and qual[qpos] < args.min_baseq:
                continue
            n_obs += 1
            base = seq[qpos]
            if base == allele[0]:
                votes[key][0] += 1
                n_hit += 1
            elif base == allele[1]:
                votes[key][1] += 1
                n_hit += 1
        if key in votes:
            umi_bin[key] = cb

    # Majority vote per UMI; a tie means the UMI saw both alleles, which is
    # either a chimera or an alignment artefact - drop it rather than guess.
    per_bin = collections.defaultdict(lambda: [0, 0])
    ties = 0
    for key, (r, a) in votes.items():
        if r > a:
            per_bin[umi_bin[key]][0] += 1
        elif a > r:
            per_bin[umi_bin[key]][1] += 1
        else:
            ties += 1

    stats["snp_observations"] += n_obs
    stats["snp_matching_an_allele"] += n_hit
    stats["umi_ties"] += ties
    return per_bin


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True)
    p.add_argument("--spatial-dir", required=True,
                   help="binned_outputs/square_002um/spatial")
    p.add_argument("--snps", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--tissue-out", default=None,
                   help="Gzipped list of every in-tissue 2um bin coordinate. "
                        "The sweep needs it for the coverage denominator: bins "
                        "with no informative UMI are absent from --out, and "
                        "counting only bins that have data would overstate "
                        "coverage badly at small tile sizes. "
                        "Default: <out>.tissue_bins.tsv.gz")
    p.add_argument("--x-chrom", default="chrX")
    p.add_argument("--autosome-chroms", default="chr19,chr18,chr17",
                   help="Control set. Small autosomes by default: the null only "
                        "needs depth comparable to chrX, and memory here scales "
                        "with UMIs on the largest chromosome processed.")
    p.add_argument("--snp-pos-col", type=int, default=2,
                   help="1-based column holding the SNP position (default 2).")
    p.add_argument("--snp-cols", default=None,
                   help="1-based REF,ALT columns. Omit to auto-detect.")
    p.add_argument("--swap-alleles", action="store_true",
                   help="Treat the second base as REF and the first as ALT. "
                        "Nothing in the bed states which strain is which; use "
                        "the autosomal ref fraction in the summary to decide.")
    p.add_argument("--snp-offset", type=int, default=0,
                   help="Added to the SNP position to reach 0-based. Use the "
                        "allele-match rate in the summary to check this.")
    p.add_argument("--min-mapq", type=int, default=255,
                   help="255 = STAR unique, which is what spaceranger emits.")
    p.add_argument("--min-baseq", type=int, default=20)
    p.add_argument("--gene-tag", default="GX")
    p.add_argument("--require-gene", action="store_true", default=True,
                   help="Only count reads assigned to a gene, matching what "
                        "Allelome.PRO2 does with its annotation bed.")
    p.add_argument("--no-require-gene", dest="require_gene", action="store_false")
    p.add_argument("--threads", type=int, default=4,
                   help="BGZF decompression threads. The counting loop is "
                        "single-threaded; this only speeds up reading.")
    p.add_argument("--drop-duplicates", action="store_true", default=True)
    p.add_argument("--keep-duplicates", dest="drop_duplicates", action="store_false")
    args = p.parse_args()

    autosomes = [c for c in args.autosome_chroms.split(",") if c]
    chroms = [args.x_chrom] + autosomes

    ref_i, alt_i, packed = sniff_snp_layout(args.snps, args.snp_cols)
    sys.stderr.write(
        "SNP layout: ref col %d, alt col %d%s, pos col %d, offset %+d%s\n"
        % (ref_i + 1, alt_i + 1,
           (" (two bases concatenated)" if packed == ""
            else " (packed, sep %r)" % packed if packed is not None else ""),
           args.snp_pos_col, args.snp_offset,
           "  [--swap-alleles]" if args.swap_alleles else "")
    )

    snps, snp_skipped = load_snps(args.snps, chroms, args.snp_pos_col - 1,
                                  ref_i, alt_i, packed, args.snp_offset,
                                  args.swap_alleles)
    for c in chroms:
        sys.stderr.write("  %s: %d SNPs\n" % (c, len(snps[c][0])))
    if snp_skipped:
        sys.stderr.write("  %d SNP lines skipped (unparseable or indel)\n" % snp_skipped)
    if not any(snps[c][0] for c in chroms):
        raise SystemExit("No SNPs loaded - check --snp-pos-col / chromosome naming.")

    pos = load_positions(args.spatial_dir)
    sys.stderr.write("Tissue-covered 2um bins: %d\n" % len(pos))

    bam = pysam.AlignmentFile(args.bam, "rb", threads=args.threads)
    bam_chroms = set(bam.references)
    missing = [c for c in chroms if c not in bam_chroms]
    if missing:
        raise SystemExit(
            "Not in the BAM header: %s\nBAM uses e.g. %s - chromosome naming "
            "differs between the BAM and the SNP bed."
            % (", ".join(missing), ", ".join(sorted(bam_chroms)[:5]))
        )

    stats = collections.Counter()
    x_counts = count_chrom(bam, args.x_chrom, *snps[args.x_chrom], args, stats)
    a_counts = collections.defaultdict(lambda: [0, 0])
    for c in autosomes:
        for bc, (r, a) in count_chrom(bam, c, *snps[c], args, stats).items():
            a_counts[bc][0] += r
            a_counts[bc][1] += a
    bam.close()

    tissue_out = args.tissue_out or (args.out + ".tissue_bins.tsv.gz")
    with gzip.open(tissue_out, "wt") as tf:
        tf.write("barcode\tarray_row\tarray_col\n")
        for bc, (r, c) in pos.items():
            tf.write("%s\t%d\t%d\n" % (bc, r, c))
    sys.stderr.write("Wrote %d tissue bin coordinates to %s\n"
                     % (len(pos), tissue_out))

    n_written = 0
    x_umi = a_umi = 0
    with open(args.out, "w") as out:
        out.write("barcode\tarray_row\tarray_col\tx_ref\tx_alt\ta_ref\ta_alt\n")
        for bc in set(x_counts) | set(a_counts):
            coord = pos.get(bc)
            if coord is None:          # off-tissue bin
                continue
            xr, xa = x_counts.get(bc, (0, 0))
            ar, aa = a_counts.get(bc, (0, 0))
            out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n"
                      % (bc, coord[0], coord[1], xr, xa, ar, aa))
            x_umi += xr + xa
            a_umi += ar + aa
            n_written += 1

    n_bins = len(pos)
    match_rate = (100.0 * stats["snp_matching_an_allele"] / stats["snp_observations"]
                  if stats["snp_observations"] else 0.0)

    sys.stderr.write("\n--- summary ---\n")
    sys.stderr.write("reads seen                 %d\n" % stats["reads_seen"])
    sys.stderr.write("reads passing filters      %d\n" % stats["reads_kept"])
    sys.stderr.write("read bases over a SNP      %d\n" % stats["snp_observations"])
    sys.stderr.write("  matching ref or alt      %.2f%%\n" % match_rate)
    sys.stderr.write("UMIs dropped as ties       %d\n" % stats["umi_ties"])
    sys.stderr.write("bins with informative UMIs %d / %d\n" % (n_written, n_bins))
    sys.stderr.write("chrX informative UMIs      %d   (lambda = %.4f per 2um bin)\n"
                     % (x_umi, x_umi / n_bins if n_bins else 0.0))
    sys.stderr.write("autosomal informative UMIs %d   (lambda = %.4f per 2um bin)\n"
                     % (a_umi, a_umi / n_bins if n_bins else 0.0))
    if a_umi:
        a_ref_frac = sum(v[0] for v in a_counts.values()) / a_umi
        sys.stderr.write("autosomal ref fraction     %.4f\n" % a_ref_frac)
        # Autosomes are biallelic, so this should sit at 0.5 nudged upward by
        # the mapping bias of the standard 10x B6 reference - reference-allele
        # reads align a little more readily than CAST ones. Landing clearly
        # below 0.5 instead means the two bases in the name field are the other
        # way round and REF is being read as CAST.
        if a_ref_frac < 0.49:
            sys.stderr.write(
                "\n*** WARNING: autosomal ref fraction is %.3f, below 0.5.\n"
                "*** Mapping to a B6 reference should bias this upward, so the\n"
                "*** allele order is probably reversed. Re-run with --swap-alleles.\n"
                % a_ref_frac)
        elif a_ref_frac > 0.60:
            sys.stderr.write(
                "\n*** NOTE: autosomal ref fraction is %.3f, a larger reference\n"
                "*** bias than the standard 10x B6 reference usually produces.\n"
                "*** Worth a look before trusting absolute ratios.\n" % a_ref_frac)

    sys.stderr.write("\nlambda x (s/2)^2 is the expected informative UMIs in an "
                     "s-um tile; the sweep does this properly.\n")

    if match_rate < 90.0 and stats["snp_observations"]:
        sys.stderr.write(
            "\n*** WARNING: only %.1f%% of bases over a SNP position carry either "
            "declared allele.\n*** That is what an off-by-one looks like. Re-run with "
            "--snp-offset -1 (or +1) and\n*** compare - the correct offset should be "
            "well above 95%%.\n" % match_rate
        )


if __name__ == "__main__":
    main()
