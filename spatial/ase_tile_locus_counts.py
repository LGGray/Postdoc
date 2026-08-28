#!/usr/bin/env python3
"""Per-tile allelic counts straight from the BAM, replacing the sinto +
Allelome.PRO2 route for the tile ratio map.

WHAT THIS REPLACES
------------------
The existing route splits the BAM into one file per 64 um tile with sinto, runs
Allelome.PRO2 on each, and reads the resulting locus_table.txt files back in
(spatial_sinto_tiles.slurm -> spatial_tile_map.slurm, ~14 h/sample). The tiles
are only ever aggregated to whole chromosomes afterwards: tile_ratio_map.R's
read_locus() sums `chr == "chrX"` and `chr %in% chr1..chr19` and throws gene
identity away, because Allelome.PRO2 was run against chr_annotation_mm39.bed -
whole-chromosome intervals, one locus_table row per chromosome.

So the BAM split buys nothing that a single pass cannot produce directly. This
script reads the BAM ONCE, assigns every read to its 64 um tile from the `sb`
tag, and accumulates per (tile, chromosome). Same numbers, one pass.

TWO DELIBERATE DIFFERENCES FROM ALLELOME.PRO2
---------------------------------------------
1. UMI collapse. Allelome.PRO2 is read-level, so PCR duplicates of one molecule
   count as independent evidence. This counts both ways in the same pass and
   writes both: a1_umi/a2_umi (one vote per UB, majority over its reads) and
   a1_reads/a2_reads (one vote per read, majority over its SNPs). --ap2-level
   picks which pair the Allelome.PRO2-format tree carries; the tidy table always
   carries both, so the comparison is free.

2. No gene requirement, and no exon requirement. Allelome.PRO2 counts over the
   annotation's intervals; here every SNP-overlapping read on the chromosome
   counts. On chr_annotation_mm39.bed those are the same thing, which is why
   this is a replication rather than a different measurement.

Everything else matches the existing counter's conventions, because it imports
them from it: read filters, the `sb` -> 2 um bin resolution, --snp-offset, the
majority vote, tie handling.

OUTPUTS  (--out-dir)
--------------------
  tile_chrom_counts.tsv          tile x chromosome, both levels. The artefact
                                 worth keeping - every downstream question is a
                                 groupby on this.
  tile_gene_counts.tsv.gz        tile x gene, both levels. Only with
                                 --annotation; this is what the Allelome.PRO2
                                 route could not give per tile at all.
  tile_locus_counts.provenance.tsv

  <ap2-dir>/<tile>/locus_table.txt   Allelome.PRO2-format adapter tree, so
                                 spatial/tile_ratio_map.R draws these counts
                                 UNMODIFIED. Written under
                                 ase/<sample>/Allelome.PRO2_tiles_<um>um_pysam
                                 because that script resolves its input
                                 directory from BASE and SNP_LABEL. The name is
                                 new, so no existing tree is touched.

WHY THE ADAPTER TREE RATHER THAN A NEW PLOTTING SCRIPT
------------------------------------------------------
tile_ratio_map.R holds the lab's allelic-ratio colour ramp, the H&E overlay, the
overdispersed autosomal null and the call logic, all tuned and all commented.
Reimplementing it would produce figures that are similar to the old ones and not
comparable to them. Emitting its input format instead means the only thing that
changed between the two sets of PDFs is where the counts came from.

USAGE
-----
  python3 ase_tile_locus_counts.py \
      --bam  <sample>/outs/possorted_genome_bam.bam \
      --spatial-dir <sample>/outs/binned_outputs/square_002um/spatial \
      --snps SNPfile_..._mm39_no_Xist.bed \
      --out-dir ase_pysam_64um/<sample> \
      --ap2-dir ase/<sample>/Allelome.PRO2_tiles_64um_pysam \
      --tile-um 64 --snp-offset -1 --min-mapq 255 --min-baseq 20 --threads 8

Normally driven by slurm/spatial_tile_locus_map.slurm, which also runs the plots.
"""

import argparse
import bisect
import collections
import datetime
import gzip
import os
import shutil
import sys

import pysam

# The counting conventions live in the sibling module and are imported rather
# than restated: a second copy of the read filters or of the SNP offset handling
# would drift from the first one silently, and the whole point of this script is
# that its numbers are comparable to that one's.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ase_bin_allele_counts as base   # noqa: E402


def decide(r, a):
    """Majority allele, or None for a tie.

    A tie means the molecule (or read) saw both alleles in equal measure, which
    is a chimera or an alignment artefact. Dropped rather than guessed, exactly
    as ase_bin_allele_counts.py does it.
    """
    return 0 if r > a else (1 if a > r else None)


def count_chrom(bam, chrom, snp_sorted, snp_map, args, stats, resolve, k,
                tiles, genes=None, gene_iv=None):
    """Accumulate one chromosome into tiles[(trow, tcol, chrom)].

    Each cell is [a1_umi, a2_umi, a1_reads, a2_reads]. Read-level counts land
    during the pass; UMI-level ones need the whole chromosome first, since a
    molecule's reads are scattered across it.
    """
    votes = {}          # (bin, umi) -> [ref, alt, leftmost informative pos]
    umi_tile = {}       # (bin, umi) -> (trow, tcol)
    tile_cache = {}     # bin barcode -> (trow, tcol) or None, off tissue

    def tile_of(bc):
        t = tile_cache.get(bc, False)
        if t is False:
            c = resolve(bc)
            t = None if c is None else (int(c[0]) // k, int(c[1]) // k)
            tile_cache[bc] = t
        return t

    for read in bam.fetch(chrom):
        stats["reads_seen"] += 1
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            stats["drop_not_primary"] += 1
            continue
        if not args.keep_duplicates and read.is_duplicate:
            stats["drop_duplicate"] += 1
            continue
        if read.mapping_quality < args.min_mapq:
            stats["drop_mapq"] += 1
            continue
        if not read.has_tag(args.bin_tag):
            stats["drop_no_bin_tag"] += 1
            continue
        bin_bc = read.get_tag(args.bin_tag)
        if not bin_bc.startswith(base.BIN_PREFIX):
            # The barcode failed correction, so `sb` holds the raw sequence and
            # the read has no position on the slide at all.
            stats["drop_bin_unresolved"] += 1
            continue
        if not read.has_tag("UB"):
            stats["drop_no_UB"] += 1
            continue
        t = tile_of(bin_bc)
        if t is None:
            stats["drop_off_tissue"] += 1
            continue
        stats["reads_kept"] += 1

        # One binary search decides whether this read can be informative at all.
        # get_aligned_pairs() builds a Python list per read and is by far the
        # most expensive call here, so it is only ever reached by reads that
        # actually span a SNP.
        lo = bisect.bisect_left(snp_sorted, read.reference_start)
        if lo >= len(snp_sorted) or snp_sorted[lo] >= read.reference_end:
            continue

        seq = read.query_sequence
        qual = read.query_qualities
        r_ref = r_alt = 0
        left = None
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            allele = snp_map.get(rpos)
            if allele is None:
                continue
            if qual is not None and qual[qpos] < args.min_baseq:
                continue
            stats["snp_observations"] += 1
            b = seq[qpos]
            if b == allele[0]:
                r_ref += 1
            elif b == allele[1]:
                r_alt += 1
            else:
                continue
            stats["snp_matching_an_allele"] += 1
            if left is None or rpos < left:
                left = rpos
        if left is None:
            continue

        cell = tiles.get((t[0], t[1], chrom))
        if cell is None:
            cell = tiles[(t[0], t[1], chrom)] = [0, 0, 0, 0]

        # ---- read level: one vote per read, majority over its own SNPs ----
        rb = decide(r_ref, r_alt)
        if rb is None:
            stats["read_ties"] += 1
        else:
            cell[2 + rb] += 1
            if genes is not None:
                nm = gene_iv.name_at(chrom, left)
                if nm is not None:
                    g = genes.get((t[0], t[1], chrom, nm))
                    if g is None:
                        g = genes[(t[0], t[1], chrom, nm)] = [0, 0, 0, 0]
                    g[2 + rb] += 1

        # ---- UMI level: accumulate this read's evidence onto the molecule ----
        key = (bin_bc, read.get_tag("UB"))
        v = votes.get(key)
        if v is None:
            v = votes[key] = [0, 0, left]
            umi_tile[key] = t
        v[0] += r_ref
        v[1] += r_alt
        # The leftmost informative SNP of the whole molecule, not of this read.
        # On 3' nuclear data a read can start a long way from the SNP that
        # carries its allele, so the SNP is what names the gene.
        if left < v[2]:
            v[2] = left

    for key, (r, a, pos) in votes.items():
        b = decide(r, a)
        if b is None:
            stats["umi_ties"] += 1
            continue
        t = umi_tile[key]
        cell = tiles.get((t[0], t[1], chrom))
        if cell is None:
            cell = tiles[(t[0], t[1], chrom)] = [0, 0, 0, 0]
        cell[b] += 1
        if genes is not None:
            nm = gene_iv.name_at(chrom, pos)
            if nm is not None:
                g = genes.get((t[0], t[1], chrom, nm))
                if g is None:
                    g = genes[(t[0], t[1], chrom, nm)] = [0, 0, 0, 0]
                g[b] += 1


def tile_name(tile_um, trow, tcol):
    """Exactly the name spatial/ase_tile_sweep.R writes into the sinto map, and
    the one tile_ratio_map.R parses its row/col back out of. Any drift here
    silently unmatches every tile."""
    return "tile_%dum_r%04d_c%04d" % (tile_um, trow, tcol)


def write_ap2_tree(out_dir, tiles, tile_um, chroms, level, force):
    """One Allelome.PRO2-format locus_table.txt per tile.

    tile_ratio_map.R's read_locus() needs columns chr, a1_reads, a2_reads and
    total_reads (it lowercases the header, so case does not matter) and requires
    a chrX row to accept the tile. A tile with no informative chrX molecule is
    therefore skipped rather than written with zeros: written, it would land in
    the map as "pending" - the colour that means "sinto has not scored this
    yet" - which would be a lie about a tile that has been scored and is empty.
    """
    if os.path.isdir(out_dir):
        if not force:
            raise SystemExit(
                "Refusing to write into an existing %s\n"
                "  Pass --force to replace it. It is a generated adapter tree, "
                "but an\n  Allelome.PRO2 output tree looks identical from here "
                "and must not be\n  overwritten by accident." % out_dir)
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    i0, i1 = (0, 1) if level == "umi" else (2, 3)
    per_tile = collections.defaultdict(dict)
    for (trow, tcol, chrom), c in tiles.items():
        per_tile[(trow, tcol)][chrom] = c

    n_written = n_skipped = 0
    for (trow, tcol), by_chrom in per_tile.items():
        x = by_chrom.get("chrX")
        if x is None or (x[i0] + x[i1]) == 0:
            n_skipped += 1
            continue
        nm = tile_name(tile_um, trow, tcol)
        d = os.path.join(out_dir, nm)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "locus_table.txt"), "w") as fh:
            fh.write("chr\tname\ta1_reads\ta2_reads\ttotal_reads\t"
                     "allelic_ratio\n")
            for chrom in chroms:
                c = by_chrom.get(chrom)
                if c is None:
                    continue
                a1, a2 = c[i0], c[i1]
                n = a1 + a2
                if n == 0:
                    continue
                fh.write("%s\t%s\t%d\t%d\t%d\t%.6f\n"
                         % (chrom, chrom, a1, a2, n, a1 / n))
        n_written += 1
    return n_written, n_skipped


def main():
    p = argparse.ArgumentParser(
        description="Per-tile allelic counts from the BAM, replacing "
                    "sinto + Allelome.PRO2 for the tile ratio map.")
    p.add_argument("--bam", required=True)
    p.add_argument("--spatial-dir", required=True,
                   help="square_002um/spatial, for tissue_positions")
    p.add_argument("--snps", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--ap2-dir", default=None,
                   help="Where to write the Allelome.PRO2-format adapter tree. "
                        "Omit to write only the tidy tables.")
    p.add_argument("--tile-um", type=int, default=64,
                   help="Must be an even number of microns: tiles are whole "
                        "2 um bins per side so the tiling nests exactly.")
    p.add_argument("--chroms", default=None,
                   help="Default chr1..chr19,chrX - chr1..chr19 because that "
                        "is the autosomal control tile_ratio_map.R pools, not "
                        "the 3-chromosome control the bin counter uses.")
    p.add_argument("--annotation", default=None,
                   help="Gene BED for the optional tile x gene table. Real BED "
                        "coordinates, 0-based half-open.")
    p.add_argument("--ap2-level", choices=("umi", "reads"), default="umi",
                   help="Which count pair the adapter tree carries. 'reads' "
                        "reproduces Allelome.PRO2's read-level statistic; "
                        "'umi' (default) is the better measurement.")
    p.add_argument("--snp-offset", type=int, default=-1)
    p.add_argument("--snp-pos-col", type=int, default=2)
    p.add_argument("--snp-cols", default=None)
    p.add_argument("--swap-alleles", action="store_true")
    p.add_argument("--min-mapq", type=int, default=255)
    p.add_argument("--min-baseq", type=int, default=20)
    p.add_argument("--keep-duplicates", action="store_true",
                   help="Off by default: PCR duplicates are dropped, matching "
                        "the bin counter.")
    p.add_argument("--bin-tag", default="sb")
    p.add_argument("--barcode-map", default=None)
    p.add_argument("--threads", type=int, default=4,
                   help="BGZF decompression threads. The counting loop itself "
                        "is single-threaded Python.")
    p.add_argument("--force", action="store_true",
                   help="Replace an existing --ap2-dir.")
    args = p.parse_args()

    if args.tile_um % 2:
        raise SystemExit("--tile-um must be even: a tile is a whole number of "
                         "2 um bins per side.")
    k = args.tile_um // 2

    chroms = (args.chroms.split(",") if args.chroms else
              ["chr%d" % i for i in range(1, 20)] + ["chrX"])
    os.makedirs(args.out_dir, exist_ok=True)

    ref_i, alt_i, packed = base.sniff_snp_layout(args.snps, args.snp_cols)
    sys.stderr.write("SNP layout: ref col %d, alt col %d%s, pos col %d, "
                     "offset %d%s\n"
                     % (ref_i + 1, alt_i + 1,
                        "" if packed is None else " (packed '%s')" % packed,
                        args.snp_pos_col, args.snp_offset,
                        "  [--swap-alleles]" if args.swap_alleles else ""))
    snps, skipped = base.load_snps(args.snps, chroms, args.snp_pos_col - 1,
                                   ref_i, alt_i, packed, args.snp_offset,
                                   args.swap_alleles)
    total_snps = sum(len(snps[c][0]) for c in chroms)
    if not total_snps:
        raise SystemExit("No SNPs loaded - check --snp-pos-col and that the "
                         "bed's chromosome names match the BAM's.")
    sys.stderr.write("SNPs loaded: %d over %d chromosomes (%d unusable lines)\n"
                     % (total_snps, len(chroms), skipped))

    gene_iv = None
    if args.annotation:
        gene_iv = base.Intervals(args.annotation)
        sys.stderr.write("Gene annotation: %s\n" % args.annotation)

    pos = base.load_positions(args.spatial_dir)
    sys.stderr.write("Tissue-covered 2 um bins: %d\n" % len(pos))
    bcmap = None
    if args.barcode_map:
        bcmap, _ = base.load_barcode_map(args.barcode_map, args.spatial_dir, pos)
    pos_stripped = {b.rsplit("-", 1)[0]: v for b, v in pos.items()}

    def resolve(bc):
        c = pos.get(bc)
        if c is not None:
            return c
        c = pos_stripped.get(bc.rsplit("-", 1)[0])
        if c is not None:
            return c
        if bcmap is not None:
            return bcmap.get(bc) or bcmap.get(bc.rsplit("-", 1)[0])
        return None

    bam = pysam.AlignmentFile(args.bam, "rb", threads=args.threads)
    bam_chroms = set(bam.references)
    missing = [c for c in chroms if c not in bam_chroms]
    if missing:
        raise SystemExit(
            "Not in the BAM header: %s\nBAM uses e.g. %s - chromosome naming "
            "differs between the BAM and the SNP bed."
            % (", ".join(missing), ", ".join(sorted(bam_chroms)[:5])))

    # Pre-flight. This pass takes hours; an unresolvable barcode namespace is
    # otherwise only discovered at the end, as an empty table.
    peek = set()
    for read in bam.fetch("chrX"):
        if read.has_tag(args.bin_tag):
            v = read.get_tag(args.bin_tag)
            if v.startswith(base.BIN_PREFIX):
                peek.add(v)
        if len(peek) >= 200:
            break
    n_ok = sum(1 for b in peek if resolve(b) is not None)
    sys.stderr.write("Pre-flight: %d / %d sampled bin tags are on tissue\n"
                     % (n_ok, len(peek)))
    if peek and n_ok == 0:
        raise SystemExit(
            "%s tags cannot be resolved to 2 um bins.\nBAM examples:\n  %s\n"
            "tissue_positions barcodes:\n  %s\nPass --barcode-map with "
            "spaceranger's sequence -> bin table."
            % (args.bin_tag, "\n  ".join(sorted(peek)[:3]),
               "\n  ".join(sorted(pos)[:3])))

    tiles = {}
    genes = {} if gene_iv is not None else None
    stats = collections.Counter()
    started = datetime.datetime.now()
    for chrom in chroms:
        t0 = datetime.datetime.now()
        snp_sorted, snp_map = snps[chrom]
        if not snp_sorted:
            sys.stderr.write("  %-6s no SNPs, skipped\n" % chrom)
            continue
        count_chrom(bam, chrom, snp_sorted, snp_map, args, stats, resolve, k,
                    tiles, genes, gene_iv)
        el = (datetime.datetime.now() - t0).total_seconds()
        sys.stderr.write("  %-6s %8d SNPs   %6.1f s   (%d tile rows so far)\n"
                         % (chrom, len(snp_sorted), el, len(tiles)))
        sys.stderr.flush()
    bam.close()
    elapsed = (datetime.datetime.now() - started).total_seconds()

    # ------------------------------------------------------------ tidy table
    tsv = os.path.join(args.out_dir, "tile_chrom_counts.tsv")
    with open(tsv, "w") as fh:
        fh.write("tile\ttrow\ttcol\tchrom\ta1_umi\ta2_umi\tn_umi\t"
                 "a1_reads\ta2_reads\tn_reads\n")
        for (trow, tcol, chrom) in sorted(tiles):
            c = tiles[(trow, tcol, chrom)]
            fh.write("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n"
                     % (tile_name(args.tile_um, trow, tcol), trow, tcol, chrom,
                        c[0], c[1], c[0] + c[1], c[2], c[3], c[2] + c[3]))
    sys.stderr.write("\nWrote %d tile x chromosome rows to %s\n"
                     % (len(tiles), tsv))

    if genes is not None:
        gz = os.path.join(args.out_dir, "tile_gene_counts.tsv.gz")
        with gzip.open(gz, "wt") as fh:
            fh.write("tile\ttrow\ttcol\tchrom\tgene\ta1_umi\ta2_umi\tn_umi\t"
                     "a1_reads\ta2_reads\tn_reads\n")
            for key in sorted(genes):
                trow, tcol, chrom, nm = key
                c = genes[key]
                fh.write("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n"
                         % (tile_name(args.tile_um, trow, tcol), trow, tcol,
                            chrom, nm, c[0], c[1], c[0] + c[1],
                            c[2], c[3], c[2] + c[3]))
        sys.stderr.write("Wrote %d tile x gene rows to %s\n" % (len(genes), gz))

    n_written = n_skipped = 0
    if args.ap2_dir:
        n_written, n_skipped = write_ap2_tree(args.ap2_dir, tiles,
                                              args.tile_um, chroms,
                                              args.ap2_level, args.force)
        sys.stderr.write("Wrote %d locus_table.txt files to %s (%s level); "
                         "%d tiles had no informative chrX molecule\n"
                         % (n_written, args.ap2_dir, args.ap2_level, n_skipped))

    # ------------------------------------------------------------- provenance
    prov = os.path.join(args.out_dir, "tile_locus_counts.provenance.tsv")
    x_umi = sum(c[0] + c[1] for (r, cc, ch), c in tiles.items() if ch == "chrX")
    x_a1 = sum(c[0] for (r, cc, ch), c in tiles.items() if ch == "chrX")
    a_umi = sum(c[0] + c[1] for (r, cc, ch), c in tiles.items() if ch != "chrX")
    a_a1 = sum(c[0] for (r, cc, ch), c in tiles.items() if ch != "chrX")
    x_rd = sum(c[2] + c[3] for (r, cc, ch), c in tiles.items() if ch == "chrX")
    x_rd1 = sum(c[2] for (r, cc, ch), c in tiles.items() if ch == "chrX")
    with open(prov, "w") as fh:
        fh.write("key\tvalue\n")
        for kk, vv in (
            ("argv", " ".join(sys.argv)),
            ("run_at", started.strftime("%Y-%m-%d %H:%M:%S")),
            ("elapsed_s", "%.0f" % elapsed),
            ("bam", args.bam),
            ("spatial_dir", args.spatial_dir),
            ("snp_bed", args.snps),
            ("snp_bed_md5", base.md5_of(args.snps)),
            ("annotation", args.annotation or "(none)"),
            ("annotation_md5",
             base.md5_of(args.annotation) if args.annotation else "(none)"),
            ("tile_um", args.tile_um),
            ("chroms", ",".join(chroms)),
            ("snp_offset", args.snp_offset),
            ("min_mapq", args.min_mapq),
            ("min_baseq", args.min_baseq),
            ("drop_duplicates", not args.keep_duplicates),
            ("ap2_dir", args.ap2_dir or "(none)"),
            ("ap2_level", args.ap2_level),
            ("tiles_with_locus_table", n_written),
            ("tiles_without_chrX", n_skipped),
            ("tile_chrom_rows", len(tiles)),
            ("tile_gene_rows", len(genes) if genes is not None else "(none)"),
            ("chrX_informative_umis", x_umi),
            ("chrX_pooled_a1_frac", "%.4f" % (x_a1 / x_umi) if x_umi else "NA"),
            ("chrX_informative_reads", x_rd),
            ("chrX_pooled_a1_frac_reads",
             "%.4f" % (x_rd1 / x_rd) if x_rd else "NA"),
            ("autosomal_informative_umis", a_umi),
            ("autosomal_pooled_a1_frac",
             "%.4f" % (a_a1 / a_umi) if a_umi else "NA"),
        ):
            fh.write("%s\t%s\n" % (kk, vv))
    sys.stderr.write("Wrote %s\n" % prov)

    # ---------------------------------------------------------------- summary
    sys.stderr.write("\n--- summary (%.1f min) ---\n" % (elapsed / 60.0))
    sys.stderr.write("reads seen                 %d\n" % stats["reads_seen"])
    for kk, lab in (("drop_not_primary", "not primary"),
                    ("drop_duplicate", "PCR duplicate"),
                    ("drop_mapq", "MAPQ below cutoff"),
                    ("drop_no_bin_tag", "no %s tag" % args.bin_tag),
                    ("drop_bin_unresolved",
                     "%s is not a 2 um bin" % args.bin_tag),
                    ("drop_no_UB", "no UB tag"),
                    ("drop_off_tissue", "bin is off tissue")):
        if stats[kk]:
            sys.stderr.write("  dropped, %-28s %d\n" % (lab, stats[kk]))
    sys.stderr.write("reads passing filters      %d\n" % stats["reads_kept"])
    sys.stderr.write("read bases over a SNP      %d  (%d matched an allele)\n"
                     % (stats["snp_observations"],
                        stats["snp_matching_an_allele"]))
    sys.stderr.write("reads dropped as ties      %d\n" % stats["read_ties"])
    sys.stderr.write("UMIs dropped as ties       %d\n" % stats["umi_ties"])
    sys.stderr.write("chrX informative UMIs      %d   pooled a1 fraction %s\n"
                     % (x_umi, "%.4f" % (x_a1 / x_umi) if x_umi else "NA"))
    sys.stderr.write("chrX informative reads     %d   pooled a1 fraction %s\n"
                     % (x_rd, "%.4f" % (x_rd1 / x_rd) if x_rd else "NA"))
    sys.stderr.write("autosomal informative UMIs %d   pooled a1 fraction %s\n"
                     % (a_umi, "%.4f" % (a_a1 / a_umi) if a_umi else "NA"))
    sys.stderr.write(
        "\na1 is the SNP bed's REF allele = B6, a2 = CAST. In this genotype B6\n"
        "carries the Xist deletion and cannot be inactivated, so CAST is the\n"
        "inactive X in every cell and the chrX a2 fraction (1 - a1) is escape.\n"
        "The autosomal a1 fraction above minus 0.5 is the B6-ward mapping bias,\n"
        "which makes escape UNDER-estimated.\n")


if __name__ == "__main__":
    main()
