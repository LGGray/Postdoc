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

# Spaceranger's 2um bin barcodes: s_002um_<row>_<col>-1
BIN_PREFIX = "s_002um_"


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


def load_barcode_map(path, spatial_dir, pos):
    """Raw CB sequence -> (array_row, array_col) for tissue-covered bins.

    The Visium HD BAM tags reads with the raw spatial barcode sequence
    (GATGACTGCGATG...), while every binned output is keyed by the bin name
    (s_002um_00596_01672-1). Nothing in the sequence encodes position;
    spaceranger writes the translation as outs/barcode_mappings.parquet.

    Resolved to coordinates as it is read, and restricted to bins that are on
    tissue: the full table is ~11M rows, and keeping it as Python strings would
    cost several GB to hold something we immediately throw most of away.

    Columns are identified by content rather than name - the names have moved
    between spaceranger releases, and the file may carry a bin column per bin
    size, of which only the 2um one is wanted.
    """
    outs = os.path.dirname(os.path.dirname(os.path.dirname(spatial_dir)))
    candidates = [path] if path else [
        os.path.join(outs, "barcode_mappings.parquet"),
        os.path.join(spatial_dir, "barcode_mappings.parquet"),
        os.path.join(outs, "spatial", "barcode_mappings.parquet"),
    ]
    src = next((c for c in candidates if c and os.path.exists(c)), None)
    if src is None:
        return None, candidates
    if not src.endswith(".parquet"):
        raise SystemExit("--barcode-map expects a .parquet file, got %s" % src)

    try:
        import pyarrow.parquet as pyparquet
    except ImportError:
        raise SystemExit("Need pyarrow to read %s" % src)

    pf = pyparquet.ParquetFile(src)
    head = next(pf.iter_batches(batch_size=8)).to_pydict()
    seq_col = bin_col = None
    row_col = col_col = None
    for name, vals in head.items():
        n = str(name).lower()
        if n in ("array_row", "row"):
            row_col = name
        elif n in ("array_col", "col", "column"):
            col_col = name
    for name, vals in head.items():
        if not vals or vals[0] is None:
            continue
        v = str(vals[0])
        stem = v.rsplit("-", 1)[0]
        if seq_col is None and len(stem) >= 14 and set(stem.upper()) <= set("ACGTN"):
            seq_col = name
        elif bin_col is None and v.startswith("s_002um_"):
            bin_col = name
    if bin_col is None and seq_col is not None and row_col and col_col:
        # A whitelist giving coordinates directly needs no bin names at all.
        mapping = {}
        n_rows = 0
        for batch in pf.iter_batches(batch_size=500000,
                                     columns=[seq_col, row_col, col_col]):
            dd = batch.to_pydict()
            for seq, r, c in zip(dd[seq_col], dd[row_col], dd[col_col]):
                n_rows += 1
                if seq is None or r is None or c is None:
                    continue
                mapping[str(seq).rsplit("-", 1)[0]] = (int(r), int(c))
        sys.stderr.write("Barcode map: %s\n  %s -> (%s, %s), %d rows\n"
                         % (src, seq_col, row_col, col_col, n_rows))
        return mapping, [src]

    if bin_col is None:
        raise SystemExit(
            "No 2um-bin column in %s\nColumns and first values:\n  %s"
            % (src, "\n  ".join("%s = %r" % (k, v[0] if v else None)
                                 for k, v in head.items()))
        )
    if seq_col is None:
        # outs/barcode_mappings.parquet relates 2um bins to coarser bins and to
        # segmented cells. It never carries the raw barcode sequence, so it
        # cannot translate a CB tag - that needs the slide whitelist.
        raise SystemExit(
            "%s has no barcode-sequence column, so it cannot translate CB "
            "tags.\nIt holds: %s\n\n"
            "This file maps bins to coarser bins and to cells, which is a "
            "different job.\nThe CB tag carries the raw spatial barcode "
            "sequence, and turning that into\na 2um bin needs the Visium HD "
            "slide whitelist that ships with spaceranger.\n"
            "Pass it with --barcode-map: either a two-column table of "
            "sequence and\ns_002um_* bin name, or sequence and array "
            "row/col.\n"
            % (src, ", ".join(head))
        )

    pos_stripped = {b.rsplit("-", 1)[0]: v for b, v in pos.items()}
    mapping = {}
    n_rows = n_dup = 0
    for batch in pf.iter_batches(batch_size=500000, columns=[seq_col, bin_col]):
        d = batch.to_pydict()
        for seq, b in zip(d[seq_col], d[bin_col]):
            n_rows += 1
            if seq is None or b is None:
                continue
            coord = pos.get(b)
            if coord is None:
                coord = pos_stripped.get(str(b).rsplit("-", 1)[0])
            if coord is None:          # bin is off tissue
                continue
            key = str(seq).rsplit("-", 1)[0]
            if key in mapping:
                n_dup += 1
                continue
            mapping[key] = coord

    sys.stderr.write("Barcode map: %s\n  %s -> %s, %d rows, %d resolve to a "
                     "tissue bin\n" % (src, seq_col, bin_col, n_rows, len(mapping)))
    if n_dup:
        # At 2um the barcodes are printed one per bin, so this should be zero.
        # If it is not, a sequence is backing more than one bin and the first
        # one wins, which is a guess - worth knowing about.
        sys.stderr.write("  WARNING: %d sequences map to more than one tissue "
                         "bin; kept the first\n" % n_dup)
    return mapping, [src]


# ------------------------------------------------------------------ main pass

def count_chrom(bam, chrom, snp_sorted, snp_map, args, stats):
    """UMI-collapsed (2um bin -> [ref, alt]) for one chromosome."""
    votes = collections.defaultdict(lambda: [0, 0])   # (bin, umi) -> votes
    umi_bin = {}                                      # (bin, umi) -> bin
    n_obs = n_hit = 0

    for read in bam.fetch(chrom):
        stats["reads_seen"] += 1
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            stats["drop_not_primary"] += 1
            continue
        if args.drop_duplicates and read.is_duplicate:
            stats["drop_duplicate"] += 1
            continue
        if read.mapping_quality < args.min_mapq:
            stats["drop_mapq"] += 1
            continue
        # Spaceranger writes the 2um bin straight onto the read as sb. When the
        # barcode failed correction it holds the raw sequence instead, and that
        # read has no position on the slide at all.
        if not read.has_tag(args.bin_tag):
            stats["drop_no_bin_tag"] += 1
            continue
        bin_bc = read.get_tag(args.bin_tag)
        if not bin_bc.startswith(BIN_PREFIX):
            stats["drop_bin_unresolved"] += 1
            continue
        if not read.has_tag("UB"):
            stats["drop_no_UB"] += 1
            continue
        if args.require_gene and not read.has_tag(args.gene_tag):
            stats["drop_no_gene"] += 1
            continue
        stats["reads_kept"] += 1

        # Most reads in a CAST cross do overlap a SNP, but skipping the ones
        # that cannot still avoids building an aligned-pairs list per read.
        lo = bisect.bisect_left(snp_sorted, read.reference_start)
        if lo >= len(snp_sorted) or snp_sorted[lo] >= read.reference_end:
            continue

        key = (bin_bc, read.get_tag("UB"))
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
            umi_bin[key] = bin_bc

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


def probe_offsets(bam, chrom, snp_sorted, snp_map, args, candidates=(-1, 0, 1)):
    """Which coordinate convention does the SNP bed use?

    A bed written from a VCF with `pos, pos+1` has a 1-based column 2, so the
    0-based position pysam reports is one less. Guessing costs a full run each
    time; this reads a slice of one chromosome and scores every candidate at
    once. The right offset puts nearly every base over a SNP on one of the two
    declared alleles - in an F1 het there is nothing else it can be, bar
    sequencing error. A wrong offset lands near the base composition, ~50%.
    """
    hit = {o: 0 for o in candidates}
    obs = {o: 0 for o in candidates}
    cb_seen = set()
    n = 0
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < args.min_mapq:
            continue
        n += 1
        if n > args.probe:
            break
        if len(cb_seen) < 5 and read.has_tag(args.bin_tag):
            cb_seen.add(read.get_tag(args.bin_tag))
        seq = read.query_sequence
        qual = read.query_qualities
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if qual is not None and qual[qpos] < args.min_baseq:
                continue
            base = seq[qpos]
            for o in candidates:
                # Stored key would be filepos + o, so the file position that
                # would line up with this base is rpos - o.
                allele = snp_map.get(rpos - o)
                if allele is None:
                    continue
                obs[o] += 1
                if base == allele[0] or base == allele[1]:
                    hit[o] += 1

    sys.stderr.write("\n--- offset probe (%d reads on %s) ---\n" % (n - 1, chrom))
    best, best_rate = None, -1.0
    for o in candidates:
        rate = 100.0 * hit[o] / obs[o] if obs[o] else 0.0
        sys.stderr.write("  --snp-offset %+d : %8d bases over a SNP, %5.1f%% "
                         "carry a declared allele\n" % (o, obs[o], rate))
        if rate > best_rate:
            best, best_rate = o, rate
    if best_rate >= 90:
        sys.stderr.write("\nUse --snp-offset %+d\n" % best)
    else:
        sys.stderr.write(
            "\nNo offset clears 90%%. The mismatch is not an off-by-one - check\n"
            "that the SNP bed and the BAM are on the same assembly (mm39 vs mm10\n"
            "would look exactly like this).\n")
    return best, sorted(cb_seen)


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
    p.add_argument("--bin-tag", default="sb",
                   help="Tag carrying the 2um bin barcode. Spaceranger writes "
                        "sb:Z:s_002um_XXXXX_YYYYY-1; CB holds the raw sequence "
                        "and cannot be placed on the slide without a whitelist.")
    p.add_argument("--gene-tag", default="GX")
    p.add_argument("--require-gene", action="store_true", default=False,
                   help="Only count reads carrying a gene tag. Off by default: "
                        "GX marks exonic assignment, whereas Allelome.PRO2 "
                        "counts within gene-body intervals from its annotation "
                        "bed, introns included. Requiring GX would be stricter "
                        "than the pipeline this is meant to be sizing, and on "
                        "a nuclear-heavy 3' assay it discards most reads.")
    p.add_argument("--barcode-map", default=None,
                   help="Spaceranger's raw-sequence -> 2um-bin table. Found "
                        "automatically in the usual places; pass it explicitly "
                        "if it lives somewhere else.")
    p.add_argument("--probe", type=int, default=0,
                   help="Read this many alignments, report the allele-match "
                        "rate for offsets -1/0/+1, and exit without counting. "
                        "200000 is plenty and takes a minute.")
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

    # The probe does its own offset arithmetic against unshifted keys, so it
    # must not be handed a pre-shifted table.
    load_offset = 0 if args.probe else args.snp_offset
    snps, snp_skipped = load_snps(args.snps, chroms, args.snp_pos_col - 1,
                                  ref_i, alt_i, packed, load_offset,
                                  args.swap_alleles)
    for c in chroms:
        sys.stderr.write("  %s: %d SNPs\n" % (c, len(snps[c][0])))
    if snp_skipped:
        sys.stderr.write("  %d SNP lines skipped (unparseable or indel)\n" % snp_skipped)
    if not any(snps[c][0] for c in chroms):
        raise SystemExit("No SNPs loaded - check --snp-pos-col / chromosome naming.")

    pos = load_positions(args.spatial_dir)
    sys.stderr.write("Tissue-covered 2um bins: %d\n" % len(pos))

    # Only consulted if --barcode-map is passed explicitly. Not auto-searched:
    # outs/barcode_mappings.parquet looks like the right file but maps bins to
    # coarser bins and to cells, and carries no barcode sequence at all.
    bcmap, tried = ((None, []) if not args.barcode_map else
                    load_barcode_map(args.barcode_map, args.spatial_dir, pos))
    pos_stripped = {b.rsplit("-", 1)[0]: v for b, v in pos.items()}

    def resolve(bc):
        """Bin barcode -> (array_row, array_col), or None if it is off tissue."""
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
            % (", ".join(missing), ", ".join(sorted(bam_chroms)[:5]))
        )

    # Pre-flight. The counting pass takes hours; resolving a handful of CB tags
    # takes a second, and an unresolvable barcode namespace is otherwise only
    # discovered at the very end, as an empty output table.
    if not args.probe:
        peek = set()
        for read in bam.fetch(args.x_chrom):
            if read.has_tag(args.bin_tag):
                v = read.get_tag(args.bin_tag)
                if v.startswith(BIN_PREFIX):
                    peek.add(v)
            if len(peek) >= 200:
                break
        n_ok = sum(1 for b in peek if resolve(b) is not None)
        sys.stderr.write("Pre-flight: %d / %d sampled bin tags are on tissue\n"
                         % (n_ok, len(peek)))
        if peek and n_ok == 0:
            sys.stderr.write(
                "\nERROR: %s tags cannot be resolved to 2um bins.\n"
                "BAM %s examples:\n  %s\nparquet barcodes:\n  %s\n"
                "%s\nPass --barcode-map with spaceranger's sequence -> bin table.\n"
                % (args.bin_tag, args.bin_tag,
                   "\n  ".join(sorted(peek)[:3]), "\n  ".join(sorted(pos)[:3]),
                   ("Barcode map loaded but no entry matched."
                    if bcmap is not None else
                    "Every sampled tag looks like a bin, so these are bins "
                    "that are\noff tissue - check the spatial dir matches "
                    "this sample."))
            )
            raise SystemExit(2)

    if args.probe:
        _, cb_seen = probe_offsets(bam, args.x_chrom, *snps[args.x_chrom], args)
        bam.close()
        # The other thing that silently produces an empty table: CB tags and
        # parquet barcodes drawn from different namespaces. Show both.
        sys.stderr.write("\n--- barcode namespace ---\n")
        sys.stderr.write("BAM %s tags:\n  %s\n"
                         % (args.bin_tag,
                            "\n  ".join(cb_seen) if cb_seen else "(none seen)"))
        sys.stderr.write("parquet barcodes:\n  %s\n"
                         % "\n  ".join(sorted(pos)[:5]))
        n_ok = sum(1 for b in cb_seen if resolve(b) is not None)
        sys.stderr.write("%d / %d sampled tags resolve to a tissue bin\n"
                         % (n_ok, len(cb_seen)))
        return

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

    # 10x writes the CB tag with a "-1" GEM-well suffix in some outputs and
    # without it in others, and the parquet does not always agree with the BAM.
    # Try the literal barcode first, then the suffix-stripped form, and refuse
    # to write an empty table quietly if neither lands.
    observed = set(x_counts) | set(a_counts)
    n_match = sum(1 for b in observed if resolve(b) is not None)
    if observed and n_match == 0:
        sys.stderr.write(
            "\nERROR: none of the %d barcodes seen in the BAM match any bin in\n"
            "%s\n\nBAM CB examples:\n  %s\n\nparquet barcode examples:\n  %s\n\n"
            "These are different barcode namespaces. If the BAM carries raw\n"
            "sequence barcodes rather than the s_002um_XXXXX_YYYYY form, the 2um\n"
            "bin identity is not recoverable from CB alone and the tiling has to\n"
            "be driven off a different tag - check what spaceranger wrote:\n"
            "  samtools view %s | head -1 | tr '\\t' '\\n' | grep -E '^(CB|UB|GX|xf)'\n"
            % (len(observed), args.spatial_dir,
               "\n  ".join(sorted(observed)[:5]),
               "\n  ".join(sorted(pos)[:5]),
               args.bam)
        )
        raise SystemExit(2)
    if observed:
        sys.stderr.write("barcodes matching a tissue bin: %d / %d (%.1f%%)\n"
                         % (n_match, len(observed),
                            100.0 * n_match / len(observed)))

    n_written = 0
    x_umi = a_umi = 0
    with open(args.out, "w") as out:
        out.write("barcode\tarray_row\tarray_col\tx_ref\tx_alt\ta_ref\ta_alt\n")
        for bc in observed:
            coord = resolve(bc)
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
    for k, lab in (("drop_not_primary", "not primary"),
                   ("drop_duplicate", "PCR duplicate"),
                   ("drop_mapq", "MAPQ below cutoff"),
                   ("drop_no_bin_tag", "no %s tag" % args.bin_tag),
                   ("drop_bin_unresolved", "%s is not a 2um bin" % args.bin_tag),
                   ("drop_no_UB", "no UB tag"),
                   ("drop_no_gene", "no %s tag" % args.gene_tag)):
        if stats[k]:
            sys.stderr.write("  dropped, %-18s %d\n" % (lab, stats[k]))
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
