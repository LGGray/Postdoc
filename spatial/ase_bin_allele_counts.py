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

Two flags move that, independently:

  --keep-duplicates   Stop honouring the duplicate flag. With the default
                      molecule unit this is not a way of counting duplicates -
                      they still collapse onto their UB - it is a way of getting
                      back the MOLECULES that deduplication removes outright,
                      which on this data is 26-38% of chrX (spaceranger flags
                      every read of some molecules as a duplicate of another
                      molecule at the same position). Those molecules are
                      independent observations and belong in the count.

  --count-unit read   Vote once per read instead of once per molecule, which is
                      Allelome.PRO2's statistic. Together with
                      --keep-duplicates it is the pre-deduplication view.
                      Use it to size tiles against THAT pipeline's read
                      thresholds, never to gain precision: the reads of one
                      molecule are one allele call repeated.

Read mode therefore also tracks how many distinct molecules its reads came from
and reports duplication_factor = informative reads / distinct (bin, UB), in the
summary and in the provenance sidecar. ase_tile_sweep.R divides by it to get an
effective n before quoting any SE. It is deliberately not reads per
non-duplicate read: deduplication still leaves ~2.5 reads per molecule here,
because the duplicate flag is per alignment position and a 3' molecule spans
several, so that ratio would under-correct by a factor of 2.5.

Three optional extra outputs, all off by default so the existing output path is
byte-for-byte unchanged:

  --subset-bed NAME=PATH        Extra count columns NAME_ref / NAME_alt, holding
                                only the UMIs whose informative SNPs fall inside
                                PATH's intervals. Repeatable. Append
                                ":complement" to the path to count the UMIs
                                OUTSIDE the intervals instead, on the
                                chromosomes the bed names.
                                A subset bed may name chromosomes that are
                                neither --x-chrom nor in --autosome-chroms; those
                                are fetched interval-by-interval, so an autosomal
                                positive control costs seconds, not a whole
                                chromosome of memory.
                                Subsets are counted in the SAME pass as the main
                                columns, so all of them partition one consistent
                                set of UMI calls - which is the point: escape vs
                                non-escape and Xic-masked vs not have to be
                                measured on the same molecules to be comparable.
                                A subset and its complement therefore sum to the
                                main column up to the handful of molecules that
                                cover informative SNPs on both sides of an
                                interval boundary, which are counted in both.

  --window-out PATH             Sparse per-window count table, allele-split:
                                chrom, win_start, tile_row, tile_col, ref, alt.
                                Windows are --window-size bp (default 100000);
                                the tile grid is --window-tile-um (default 64,
                                use 2 for one row per 2um bin). This is what
                                answers "where on the chromosome are the excess
                                chrX UMIs" - the main table has already summed
                                over position and cannot.

  --gene-bin-out PATH           Per-gene, per-pixel allelic counts on a
                                --gene-bin-um grid: chrom, gene, row, col, ref,
                                alt. This is the one output that keeps gene
                                identity AND position, which is what spASE's
                                scase()/spase() need as matrix1/matrix2
                                (spatial/spase_scase.R). Pair it with
                                --gene-bed; molecules that fall in no gene body
                                are absent here but still counted in the main
                                table.

  --gene-bed PATH               Gene-body bed (gene name in column 4) that
                                assigns molecules to genes for --gene-bin-out.
                                Use it. The alternative is the GX tag, which
                                spaceranger sets only on exonic assignment,
                                and this data is mostly intronic - the same
                                reason spatial_ase_sweep.slurm does not pass
                                --require-gene. With a bed, the gene is the
                                interval containing the molecule's leftmost
                                informative SNP, matching both the locus table
                                here and Allelome.PRO2's gene-body scoring.

  --gene-bin-um N               Pixel size for --gene-bin-out, in microns.
                                Default 16, i.e. 8x8 capture bins. Must be even.
                                spASE fits a beta-binomial per gene over pixels,
                                so it wants MANY pixels rather than deep ones;
                                64 would match the tile analysis but leaves only
                                a few hundred pixels per section.

  --locus-out PATH              Per-interval UMI counts for every named interval
                                in every --subset-bed, pooled over the section:
                                subset, locus, chrom, start, end, ref, alt.
                                Use it to check a positive control has usable
                                depth before reading anything into its curve.

Every run also writes <out>.provenance.tsv: the SNP bed path and its md5, the
same for each subset and interval bed, the argument vector, and the number of
distinct informative SNPs and genes per chromosome set. Three "no-Xist" SNP beds
are named in this repo and only one is live, so a result without its bed's
fingerprint next to it cannot be attributed later.

Run: see slurm/spatial_ase_sweep.slurm
"""

import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import os
import re
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


def load_snps(path, chroms, pos_col, ref_i, alt_i, packed_sep, offset, swap,
              subsets=()):
    """chrom -> (sorted position array, {pos: (ref, alt, members)}), 0-based.

    `members` is a tuple of subset indices this SNP position belongs to, empty
    for a SNP in no subset. Computed once here rather than per read: the
    counting loop touches every SNP hundreds of times and an interval search
    there would dominate the run. The tuples are interned, so N SNPs sharing a
    membership pattern cost one tuple, not N.
    """
    want = set(chroms)
    tables = {c: {} for c in chroms}
    skipped = 0
    intern = {}

    def members_for(chrom, pos):
        m = []
        for i, s in enumerate(subsets):
            if chrom not in s.chroms:
                continue
            if s.contains(chrom, pos) != s.complement:
                m.append(i)
        m = tuple(m)
        return intern.setdefault(m, m)

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
            tables[chrom][pos] = (ref, alt, members_for(chrom, pos))
    out = {}
    for c, d in tables.items():
        out[c] = (sorted(d), d)
    return out, skipped


# ----------------------------------------------------------- interval subsets

class Intervals:
    """Sorted per-chromosome intervals with names, for containment and lookup.

    COORDINATES. These beds are real BED: column 2 is a 0-based start and
    column 3 an exclusive end, so an interval covers [start, end). The SNP bed
    is NOT - column 2 there is a 1-based VCF position written as `pos, pos+1`,
    which is why the counter takes --snp-offset -1 to reach the 0-based
    coordinate pysam reports. Membership is therefore tested with the SNP
    position ALREADY shifted to 0-based against a half-open BED interval, which
    is the one combination that needs no further correction. If an interval bed
    ever arrives in the SNP bed's convention instead, pass
    --subset-bed-offset -1 rather than adjusting the intervals by hand.
    """

    def __init__(self, path, offset=0):
        self.path = path
        per = collections.defaultdict(list)
        with open_maybe_gz(path) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                try:
                    start, end = int(f[1]) + offset, int(f[2]) + offset
                except ValueError:
                    continue
                if end <= start:
                    continue
                name = f[3] if len(f) > 3 and f[3] not in ("", ".") else \
                    "%s:%d-%d" % (f[0], start, end)
                per[f[0]].append((start, end, name))
        self.by_chrom = {}
        for c, ivs in per.items():
            ivs.sort()
            # Prefix maximum of the ends, so a containment test can stop as soon
            # as no earlier interval can still reach the query. Overlapping
            # gene-body intervals are the norm in a RefSeq bed, so the cheaper
            # "previous interval only" test would be wrong here.
            starts, ends, names, maxend = [], [], [], []
            run = -1
            for s, e, n in ivs:
                starts.append(s)
                ends.append(e)
                names.append(n)
                run = max(run, e)
                maxend.append(run)
            self.by_chrom[c] = (starts, ends, names, maxend)
        self.chroms = set(self.by_chrom)
        self.n = sum(len(v[0]) for v in self.by_chrom.values())

    def contains(self, chrom, pos):
        t = self.by_chrom.get(chrom)
        if t is None:
            return False
        starts, ends, _, maxend = t
        i = bisect.bisect_right(starts, pos) - 1
        while i >= 0 and maxend[i] > pos:
            if ends[i] > pos:
                return True
            i -= 1
        return False

    def name_at(self, chrom, pos):
        """First interval (by start) containing pos, or None.

        Ties are resolved by start position, so a SNP inside two overlapping
        transcripts of the same gene is attributed once and deterministically.
        """
        t = self.by_chrom.get(chrom)
        if t is None:
            return None
        starts, ends, names, maxend = t
        i = bisect.bisect_right(starts, pos) - 1
        hit = None
        while i >= 0 and maxend[i] > pos:
            if ends[i] > pos:
                hit = names[i]
            i -= 1
        return hit

    def merged(self, chrom):
        t = self.by_chrom.get(chrom)
        if t is None:
            return []
        out = []
        for s, e in zip(t[0], t[1]):
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        return [(s, e) for s, e in out]


class Subset:
    """One extra pair of count columns, defined by a bed of intervals."""

    def __init__(self, name, path, complement, offset=0):
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
            raise SystemExit(
                "Bad --subset-bed name %r: use letters, digits and underscores, "
                "starting with a letter (it becomes a column name)." % name)
        if name in ("x", "a", "barcode", "array_row", "array_col"):
            raise SystemExit("--subset-bed name %r collides with an existing "
                             "column" % name)
        self.name = name
        self.complement = complement
        self.iv = Intervals(path, offset)
        self.path = path

    @property
    def chroms(self):
        return self.iv.chroms

    def contains(self, chrom, pos):
        return self.iv.contains(chrom, pos)


def parse_subset_args(specs, offset):
    """NAME=PATH, optionally NAME=PATH:complement."""
    out = []
    for spec in specs or ():
        if "=" not in spec:
            raise SystemExit("--subset-bed wants NAME=PATH, got %r" % spec)
        name, rest = spec.split("=", 1)
        complement = False
        if rest.endswith(":complement"):
            complement = True
            rest = rest[: -len(":complement")]
        if not os.path.exists(rest):
            raise SystemExit("--subset-bed %s: no such file %s" % (name, rest))
        out.append(Subset(name.strip(), rest, complement, offset))
    names = [s.name for s in out]
    if len(set(names)) != len(names):
        raise SystemExit("--subset-bed names must be unique: %s" % names)
    return out


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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

def count_chrom(bam, chrom, snp_sorted, snp_map, args, stats,
                subsets=(), sub_counts=None, primary=True, regions=None,
                windows=None, resolve=None, win_k=None, loci=None, seen=None,
                mol_seen=None, gene_bins=None, gene_k=None, gene_iv=None,
                snp_counts=None):
    """One chromosome, collapsed to (2um bin -> [ref, alt]).

    The unit is the molecule by default and the read under --count-unit read.
    Both go through the same vote-and-decide code below; the only difference is
    the key votes are held under and when they are flushed.

    Side outputs, all optional and all derived from the SAME per-UMI allele
    calls as the return value, so nothing here can disagree with the main table:

      sub_counts  list of dicts, one per subset, bin -> [ref, alt]. Accumulated
                  across chromosomes by the caller, since a subset bed may span
                  several.
      windows     (chrom, window index, tile row, tile col) -> [ref, alt].
      loci        (subset name, locus name) -> [ref, alt].
      gene_bins   (chrom, gene, tile row, tile col) -> [ref, alt], on the
                  --gene-bin-um grid. This is the gene x pixel pair of matrices
                  spASE takes as matrix1/matrix2 (spatial/spase_scase.R); every
                  other output here has already summed over gene identity and
                  cannot be turned back into one.
      snp_counts  (chrom, 0-based pos) -> [obs_ref, obs_alt, mol_ref, mol_alt].
                  The allele ledger, one row per SNP that was ever informative.
                  obs_* count BASE OBSERVATIONS, one per read that put a
                  declared allele over the SNP; mol_* count MOLECULES whose
                  leftmost informative SNP was this one, which is the same
                  attribution rule the window and gene tables use. Both are
                  kept because they answer different questions: obs_* is the
                  direct evidence that a SNP's two alleles are declared the
                  right way round, and needs no attribution at all, while
                  mol_* says how much of a gene's or window's signal that SNP
                  is responsible for. A SNP whose obs_alt/obs is ~1.0 on an
                  AUTOSOME cannot be biology - autosomes are biallelic here -
                  so the autosomal arm of this table is the calibration for
                  reading the chrX arm.

      gene_iv     Intervals of gene bodies. When given, a molecule's gene is
                  the interval containing its leftmost informative SNP, which
                  is the SAME rule emit_sub uses for the locus table and the
                  same gene-body definition Allelome.PRO2 scores against -
                  introns included. Without it the gene comes from the GX tag,
                  which spaceranger sets on EXONIC assignment only; on this
                  nuclear-heavy 3' data that is a minority of reads, which is
                  why spatial_ase_sweep.slurm does not pass --require-gene.
      seen        {"snps": set, "genes": set} for the provenance sidecar.
      mol_seen    read mode only: the (bin, UB) of every informative read, so
                  the duplication factor can be reads per MOLECULE. Reads per
                  non-duplicate read is not the same number and is not the one
                  the sweep needs - deduplication leaves ~2.5 reads per molecule
                  on this data, since the duplicate flag is per alignment
                  position and a molecule spans several.

    primary=False counts nothing into the returned table. That is the mode for a
    chromosome that appears only in a subset bed: its reads exist to fill the
    subset columns and must not join the autosomal control set.
    regions limits the fetch to a list of (start, end), which is what makes such
    a chromosome cheap.
    """
    votes = {}                                        # (bin, unit) -> [r, a, pos]
    umi_bin = {}                                      # (bin, unit) -> bin
    svotes = [{} for _ in subsets]
    n_obs = n_hit = 0
    # Read level uses the identical machinery with a per-read key, and flushes
    # the entry the moment the read is done. Nothing accumulates across the
    # chromosome, which is what makes a duplicate-inclusive pass fit: it has
    # several times as many reads as there are molecules, and holding a vote
    # for each of them would not.
    read_level = args.count_unit == "read"
    n_units = 0
    touched = []                                      # subsets this read hit
    snps_seen = seen["snps"] if seen else None
    genes_seen = seen["genes"] if seen else None

    if regions is None:
        it = bam.fetch(chrom)
    else:
        # A generator rather than a concatenated list: the intervals of a
        # positive-control bed are small, but a read overlapping two of them
        # would otherwise be counted twice, so the UMI keys deduplicate it.
        def _iter():
            for s, e in regions:
                for r in bam.fetch(chrom, max(0, s), e):
                    yield r
        it = _iter()

    # Vote -> table. Defined ahead of the loop because read mode calls them per
    # read; molecule mode still calls them once per UMI after the fetch. One
    # implementation, so the two units cannot drift apart.
    def decide(r, a):
        # Majority vote. A tie means the unit saw both alleles, which is either
        # a chimera or an alignment artefact - drop it rather than guess.
        return 0 if r > a else (1 if a > r else None)

    per_bin = collections.defaultdict(lambda: [0, 0])
    counters = {"ties": 0}
    tile_of = {}
    # A second bin -> tile cache, because the gene table's grid is a different
    # size from the window table's and one cache cannot serve both.
    gtile_of = {}
    umi_gene = {}                                     # (bin, unit) -> gene

    def emit_main(key, r, a, pos):
        if not primary:
            return
        b = decide(r, a)
        if b is None:
            counters["ties"] += 1
            return
        bc = umi_bin[key]
        per_bin[bc][b] += 1
        if snp_counts is not None:
            # Molecule level, attributed to the leftmost informative SNP -
            # the same `pos` the window and gene tables key on, so the three
            # tables partition the same molecules the same way and a SNP's
            # mol_* can be summed back to its window's or gene's totals.
            sc = snp_counts.get(pos)
            if sc is None:
                sc = snp_counts[pos] = [0, 0, 0, 0]
            sc[2 + b] += 1
        if gene_bins is not None:
            # pos is the molecule's leftmost informative SNP, already 0-based,
            # which is the coordinate convention Intervals expects (see its
            # docstring) - so this needs no offset of its own.
            g = gene_iv.name_at(chrom, pos) if gene_iv is not None \
                else umi_gene.get(key)
            # No gene means the molecule was intergenic, or ambiguously tagged
            # in the GX fallback. It stays in the main table - the
            # chromosome-wide aggregate is deliberately not gene-restricted -
            # and is simply absent here.
            if g is not None:
                t = gtile_of.get(bc, False)
                if t is False:
                    coord = resolve(bc)
                    t = None if coord is None else (int(coord[0]) // gene_k,
                                                    int(coord[1]) // gene_k)
                    gtile_of[bc] = t
                if t is not None:
                    gk = (chrom, g, t[0], t[1])
                    cell = gene_bins.get(gk)
                    if cell is None:
                        cell = gene_bins[gk] = [0, 0]
                    cell[b] += 1
        if windows is None:
            return
        t = tile_of.get(bc, False)
        if t is False:
            coord = resolve(bc)
            t = None if coord is None else (int(coord[0]) // win_k,
                                            int(coord[1]) // win_k)
            tile_of[bc] = t
        if t is not None:
            wk = (chrom, pos // args.window_size, t[0], t[1])
            w = windows.get(wk)
            if w is None:
                w = windows[wk] = [0, 0]
            w[b] += 1

    def emit_sub(si, key, r, a, spos):
        b = decide(r, a)
        if b is None:
            return
        s = subsets[si]
        tgt = sub_counts[si] if sub_counts is not None else None
        if tgt is not None:
            cell = tgt.get(umi_bin[key])
            if cell is None:
                cell = tgt[umi_bin[key]] = [0, 0]
            cell[b] += 1
        # Per-locus depth. Meaningless for a complement subset - "not in any of
        # these intervals" has no interval to attribute to - so those are
        # skipped rather than filed under a made-up name.
        if loci is not None and not s.complement:
            nm = s.iv.name_at(chrom, spos)
            if nm is not None:
                lk = (s.name, nm)
                cell = loci.get(lk)
                if cell is None:
                    cell = loci[lk] = [0, 0]
                cell[b] += 1

    for read in it:
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

        # The gene for the gene x pixel table, read once per read rather than
        # once per SNP. A read tagged with several genes is ambiguous and gets
        # none: such a molecule is still counted in the main table and is only
        # left out of the gene table, which is the conservative direction.
        gx = None
        if gene_bins is not None and gene_iv is None and \
                read.has_tag(args.gene_tag):
            gx = read.get_tag(args.gene_tag)
            if ";" in gx:
                gx = None

        # Most reads in a CAST cross do overlap a SNP, but skipping the ones
        # that cannot still avoids building an aligned-pairs list per read.
        lo = bisect.bisect_left(snp_sorted, read.reference_start)
        if lo >= len(snp_sorted) or snp_sorted[lo] >= read.reference_end:
            continue

        # In read mode the unit is this alignment, so the key must not collide
        # with the other reads of its molecule - a counter, not the UB tag.
        # UB is still required above, because a read without one is not a
        # counted molecule in any mode and would not be in the BAM the rest of
        # the pipeline sees either.
        if read_level:
            n_units += 1
            key = (bin_bc, n_units)
        else:
            key = (bin_bc, read.get_tag("UB"))
        seq = read.query_sequence
        qual = read.query_qualities
        v = None
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            allele = snp_map.get(rpos)
            if allele is None:
                continue
            if qual is not None and qual[qpos] < args.min_baseq:
                continue
            n_obs += 1
            base = seq[qpos]
            if base == allele[0]:
                bucket = 0
            elif base == allele[1]:
                bucket = 1
            else:
                continue
            n_hit += 1
            if snps_seen is not None:
                snps_seen.add(rpos)
            if snp_counts is not None:
                # Observation level, before any voting. A molecule with three
                # reads over this SNP contributes three, which is what makes
                # this the raw evidence for how the SNP itself behaves rather
                # than for how the molecules carrying it were resolved.
                sc = snp_counts.get(rpos)
                if sc is None:
                    sc = snp_counts[rpos] = [0, 0, 0, 0]
                sc[bucket] += 1
            if v is None:
                v = votes.get(key)
                if v is None:
                    v = votes[key] = [0, 0, rpos]
                umi_bin[key] = bin_bc
                # First informative read wins. Holding a vote per (molecule,
                # gene) would cost more memory than the ambiguity is worth;
                # 3' molecules assigned to two genes are rare.
                if gx is not None:
                    umi_gene.setdefault(key, gx)
            v[bucket] += 1
            # The UMI's position on the chromosome, for the window table. The
            # leftmost informative SNP, not the read start: it is the SNP that
            # carries the allele, and a 3' read can start a long way from it.
            if rpos < v[2]:
                v[2] = rpos
            for si in allele[2]:
                sv = svotes[si].get(key)
                if sv is None:
                    sv = svotes[si][key] = [0, 0, rpos]
                    if read_level:
                        touched.append(si)
                sv[bucket] += 1
                # The subset's own leftmost informative SNP, which is not
                # necessarily the UMI's: a molecule can cover a SNP inside an
                # escape gene and another outside it, and the per-locus table
                # has to be keyed on the one that put it in the subset.
                if rpos < sv[2]:
                    sv[2] = rpos
        if v is not None and genes_seen is not None and read.has_tag(args.gene_tag):
            genes_seen.add(read.get_tag(args.gene_tag))
        if v is not None:
            # Informative reads, how many survive the duplicate flag, and -
            # in read mode - which molecules they came from. reads / molecules
            # is the duplication factor the sweep deflates n by; the
            # nondup count is only a diagnostic, and is the wrong divisor,
            # because dedup leaves several reads per molecule.
            stats["informative_reads"] += 1
            if not read.is_duplicate:
                stats["informative_reads_nondup"] += 1
            if primary:
                stats["informative_reads_primary"] += 1
                if mol_seen is not None:
                    mol_seen.add((bin_bc, read.get_tag("UB")))
            if read_level:
                # Subsets first: emit_main deletes the bin this key maps to.
                for si in touched:
                    sv = svotes[si].pop(key)
                    emit_sub(si, key, sv[0], sv[1], sv[2])
                del touched[:]
                emit_main(key, v[0], v[1], v[2])
                del votes[key]
                del umi_bin[key]
                umi_gene.pop(key, None)

    # Molecule mode flushes here and not earlier: a UMI's reads are scattered
    # along the chromosome, so its vote is only complete once the fetch is done.
    # Read mode has already emitted and deleted every entry, so these loops run
    # over empty dicts.
    for key, (r, a, pos) in votes.items():
        emit_main(key, r, a, pos)
    for si, sub in enumerate(subsets):
        if chrom not in sub.chroms:
            continue
        for key, (r, a, spos) in svotes[si].items():
            emit_sub(si, key, r, a, spos)

    stats["snp_observations"] += n_obs
    stats["snp_matching_an_allele"] += n_hit
    stats["unit_ties"] += counters["ties"]
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
    p.add_argument("--count-unit", choices=("molecule", "read"),
                   default="molecule",
                   help="What one count is. 'molecule' (default) collapses "
                        "reads sharing (bin, UB) and votes once per molecule. "
                        "'read' votes once per read, which is Allelome.PRO2's "
                        "statistic; combined with --keep-duplicates it is the "
                        "pre-deduplication view. Read-level counts are NOT "
                        "independent observations - see the duplication factor "
                        "in the summary.")
    p.add_argument("--subset-bed", action="append", default=None,
                   metavar="NAME=PATH[:complement]",
                   help="Extra count columns NAME_ref/NAME_alt holding only the "
                        "UMIs whose informative SNPs fall inside PATH's "
                        "intervals; ':complement' counts the ones outside "
                        "instead. Repeatable. Chromosomes named only here are "
                        "fetched interval by interval and do not join the "
                        "autosomal control.")
    p.add_argument("--subset-bed-offset", type=int, default=0,
                   help="Added to subset interval starts and ends to reach "
                        "0-based half-open. 0 for a real BED, which is what the "
                        "lab's annotation beds are; -1 only if an interval bed "
                        "arrives in the SNP bed's 1-based convention.")
    p.add_argument("--window-out", default=None,
                   help="Sparse per-window allele counts: chrom, win_start, "
                        "tile_row, tile_col, ref, alt. Answers where on the "
                        "chromosome the UMIs are, which the main table cannot.")
    p.add_argument("--window-size", type=int, default=100000,
                   help="Window width in bp for --window-out (default 100000).")
    p.add_argument("--window-tile-um", type=int, default=64,
                   help="Tile grid for --window-out, in microns. 2 gives one "
                        "row per 2um bin, which is exact but large; 64 matches "
                        "the tile size the downstream analysis uses.")
    p.add_argument("--window-chroms", default=None,
                   help="Restrict --window-out to these chromosomes (default: "
                        "the X chromosome only, which is where the question is; "
                        "'all' includes the autosomal control).")
    p.add_argument("--gene-bin-out", default=None,
                   help="Per-gene per-pixel allelic counts, the gene x pixel "
                        "matrix spASE takes. See --gene-bin-um.")
    p.add_argument("--gene-bin-um", type=int, default=16,
                   help="Pixel size in microns for --gene-bin-out (default 16).")
    p.add_argument("--gene-bed", default=None,
                   help="Gene-body bed (name in column 4) used to assign each "
                        "molecule to a gene for --gene-bin-out. STRONGLY "
                        "preferred over the GX tag fallback: GX is exonic-only "
                        "and this data is mostly intronic.")
    p.add_argument("--snp-out", default=None,
                   help="Per-SNP allele ledger (gzipped tsv): chrom, pos, the "
                        "two declared alleles, base observations per allele, "
                        "molecules attributed to this SNP per allele, and the "
                        "gene body containing it when --gene-bed is given. "
                        "This is the table that says WHICH SNPs a skewed "
                        "region's signal comes from; every other output has "
                        "already summed over SNP identity.")
    p.add_argument("--snp-out-chroms", default="all",
                   help="Which chromosomes the --snp-out ledger covers: 'all' "
                        "(default, every chromosome in --x-chrom and "
                        "--autosome-chroms) or a comma list. Keep the "
                        "autosomes in: they are biallelic by construction, so "
                        "they are the only calibration for what an "
                        "artefactual per-SNP skew looks like on this data.")
    p.add_argument("--locus-out", default=None,
                   help="Per-interval UMI counts for every named interval in "
                        "every --subset-bed, pooled over the section. Check a "
                        "positive control has depth before reading its curve.")
    p.add_argument("--provenance-out", default=None,
                   help="Bed paths and md5s, the argument vector, and the "
                        "informative SNP and gene counts. Default: "
                        "<out>.provenance.tsv")
    args = p.parse_args()

    if args.window_tile_um % 2 != 0 or args.window_tile_um < 2:
        raise SystemExit("--window-tile-um must be an even number of microns "
                         "(the capture grid is 2um), got %d" % args.window_tile_um)
    if args.gene_bin_um % 2 != 0 or args.gene_bin_um < 2:
        raise SystemExit("--gene-bin-um must be an even number of microns "
                         "(the capture grid is 2um), got %d" % args.gene_bin_um)
    # --gene-bed feeds two outputs now: it assigns molecules to genes for
    # --gene-bin-out, and it names the containing gene for each row of
    # --snp-out. Either one is reason enough to load it.
    if args.gene_bed and not (args.gene_bin_out or args.snp_out):
        raise SystemExit("--gene-bed does nothing without --gene-bin-out "
                         "or --snp-out")
    if args.snp_out and not args.gene_bed:
        sys.stderr.write(
            "WARNING: --snp-out without --gene-bed leaves the ledger's `gene` "
            "column empty. The scan still works on position alone, but it "
            "cannot tell you which gene a flipped SNP sits in.\n")
    if args.gene_bin_out and not args.gene_bed:
        sys.stderr.write(
            "WARNING: --gene-bin-out without --gene-bed falls back to the %s "
            "tag, which spaceranger sets on EXONIC assignment only. On this "
            "nuclear-heavy 3' data that discards most molecules. Pass "
            "--gene-bed with a gene-body bed instead.\n" % args.gene_tag)

    autosomes = [c for c in args.autosome_chroms.split(",") if c]
    chroms = [args.x_chrom] + autosomes

    subsets = parse_subset_args(args.subset_bed, args.subset_bed_offset)
    # A subset bed may name chromosomes the main pass does not process. Those
    # get their own region-limited pass: the point of an autosomal positive
    # control is that it costs a handful of intervals, not chr7 in memory.
    sub_only = []
    for s in subsets:
        for c in sorted(s.chroms):
            if c not in chroms and c not in sub_only:
                if s.complement:
                    raise SystemExit(
                        "--subset-bed %s is a complement over %s, which is not "
                        "in --x-chrom or --autosome-chroms.\nA complement there "
                        "would mean the whole chromosome, which is not what the "
                        "flag is for - add %s to --autosome-chroms if that is "
                        "really the intent." % (s.name, c, c))
                sub_only.append(c)
    all_chroms = chroms + sub_only
    if subsets:
        for s in subsets:
            sys.stderr.write(
                "Subset %-12s %d intervals on %s%s\n  %s\n"
                % (s.name, s.iv.n, ",".join(sorted(s.chroms)),
                   "  [COMPLEMENT: SNPs outside them]" if s.complement else "",
                   s.path))
        if sub_only:
            sys.stderr.write("  region-limited passes for: %s\n"
                             % ", ".join(sub_only))

    snp_out_chroms = None
    if args.snp_out:
        if args.snp_out_chroms in (None, "", "all"):
            snp_out_chroms = set(chroms)
        else:
            snp_out_chroms = {c for c in args.snp_out_chroms.split(",") if c}
        if not snp_out_chroms & set(autosomes):
            sys.stderr.write(
                "WARNING: --snp-out covers no autosome. The ledger's whole "
                "point is that autosomal SNPs are biallelic by construction "
                "and so calibrate what a skewed chrX SNP means; without them "
                "there is nothing to compare a chrX skew against.\n")
    win_chroms = None
    if args.window_out:
        if args.window_chroms in (None, ""):
            win_chroms = {args.x_chrom}
        elif args.window_chroms == "all":
            win_chroms = set(chroms)
        else:
            win_chroms = {c for c in args.window_chroms.split(",") if c}

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
    snps, snp_skipped = load_snps(args.snps, all_chroms, args.snp_pos_col - 1,
                                  ref_i, alt_i, packed, load_offset,
                                  args.swap_alleles, subsets)
    for c in all_chroms:
        n_sub = ["%s %d" % (s.name, sum(1 for p in snps[c][0]
                                        if i in snps[c][1][p][2]))
                 for i, s in enumerate(subsets) if c in s.chroms]
        sys.stderr.write("  %s: %d SNPs%s\n"
                         % (c, len(snps[c][0]),
                            ("  (" + ", ".join(n_sub) + ")") if n_sub else ""))
    if snp_skipped:
        sys.stderr.write("  %d SNP lines skipped (unparseable or indel)\n" % snp_skipped)
    if not any(snps[c][0] for c in all_chroms):
        raise SystemExit("No SNPs loaded - check --snp-pos-col / chromosome naming.")
    for i, s in enumerate(subsets):
        n = sum(1 for c in s.chroms if c in snps
                for p in snps[c][0] if i in snps[c][1][p][2])
        if n == 0:
            sys.stderr.write(
                "\n*** WARNING: subset %s matches no SNP at all.\n"
                "*** Its columns will be all-zero, which downstream looks like "
                "a biological\n*** absence rather than a coordinate mismatch. "
                "Check the bed's assembly and\n*** that --subset-bed-offset "
                "matches its convention.\n" % s.name)

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
    missing = [c for c in all_chroms if c not in bam_chroms]
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
    sub_counts = [{} for _ in subsets]
    windows = {} if args.window_out else None
    loci = {} if args.locus_out else None
    seen = {}          # chromosome set -> informative SNP and gene identities
    # Read mode only. About a million entries on this data - the same order as
    # the vote dict molecule mode already holds - and it is the only way to say
    # how many independent molecules the read counts stand for.
    mol_seen = set() if args.count_unit == "read" else None
    win_k = args.window_tile_um // 2
    gene_bins = {} if args.gene_bin_out else None
    gene_k = args.gene_bin_um // 2
    gene_iv = None
    if args.gene_bed:
        # Same offset convention as the subset beds: real BED, half-open, no
        # shift, tested against an already-0-based SNP position.
        gene_iv = Intervals(args.gene_bed, args.subset_bed_offset)
        sys.stderr.write("Gene bodies: %d chromosomes, %d intervals from %s\n"
                         % (len(gene_iv.by_chrom),
                            sum(len(v[0]) for v in gene_iv.by_chrom.values()),
                            args.gene_bed))

    # Keyed by position only, one dict per chromosome, mirroring how snps[c] is
    # already passed in: count_chrom never sees more than one chromosome, so a
    # composite key would cost memory on every SNP to store what the caller
    # already knows.
    snp_counts = {} if args.snp_out else None

    def kw(c, primary=True, regions=None, label=None):
        lab = label or c
        seen.setdefault(lab, {"snps": set(), "genes": set()})
        return dict(subsets=subsets, sub_counts=sub_counts, primary=primary,
                    regions=regions, mol_seen=mol_seen,
                    windows=windows if (win_chroms and c in win_chroms) else None,
                    resolve=resolve, win_k=win_k, loci=loci, seen=seen[lab],
                    gene_bins=gene_bins, gene_k=gene_k, gene_iv=gene_iv,
                    snp_counts=(snp_counts.setdefault(c, {})
                                if (snp_counts is not None and primary
                                    and c in snp_out_chroms) else None))

    x_counts = count_chrom(bam, args.x_chrom, *snps[args.x_chrom], args, stats,
                           **kw(args.x_chrom))
    a_counts = collections.defaultdict(lambda: [0, 0])
    for c in autosomes:
        for bc, (r, a) in count_chrom(bam, c, *snps[c], args, stats,
                                      **kw(c, label="autosome")).items():
            a_counts[bc][0] += r
            a_counts[bc][1] += a
    # Chromosomes that exist only to fill a subset column: fetched interval by
    # interval, and counted into no main column.
    for c in sub_only:
        regions = []
        for s in subsets:
            regions.extend(s.iv.merged(c))
        regions.sort()
        merged = []
        for st, en in regions:
            if merged and st <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], en)
            else:
                merged.append([st, en])
        sys.stderr.write("  %s: %d subset regions, %.2f Mb\n"
                         % (c, len(merged),
                            sum(e - s for s, e in merged) / 1e6))
        count_chrom(bam, c, *snps[c], args, stats,
                    **kw(c, primary=False, regions=merged, label="subset_only"))
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
    for sc in sub_counts:
        observed |= set(sc)
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
    sub_umi = [0] * len(subsets)
    sub_ref = [0] * len(subsets)
    with open(args.out, "w") as out:
        # Subset columns are appended, never interleaved, so a reader that only
        # knows the original seven columns still reads this file correctly.
        head = ["barcode", "array_row", "array_col",
                "x_ref", "x_alt", "a_ref", "a_alt"]
        for s in subsets:
            head += [s.name + "_ref", s.name + "_alt"]
        out.write("\t".join(head) + "\n")
        for bc in observed:
            coord = resolve(bc)
            if coord is None:          # off-tissue bin
                continue
            xr, xa = x_counts.get(bc, (0, 0))
            ar, aa = a_counts.get(bc, (0, 0))
            row = [bc, coord[0], coord[1], xr, xa, ar, aa]
            for i, sc in enumerate(sub_counts):
                sr, sa = sc.get(bc, (0, 0))
                row += [sr, sa]
                sub_umi[i] += sr + sa
                sub_ref[i] += sr
            out.write("\t".join(str(v) for v in row) + "\n")
            x_umi += xr + xa
            a_umi += ar + aa
            n_written += 1

    if windows is not None:
        with gzip.open(args.window_out, "wt") as wf:
            wf.write("chrom\twin_start\ttile_row\ttile_col\tref\talt\n")
            for (c, wi, tr, tc), (r, a) in sorted(windows.items()):
                wf.write("%s\t%d\t%d\t%d\t%d\t%d\n"
                         % (c, wi * args.window_size, tr, tc, r, a))
        sys.stderr.write("Wrote %d window x tile rows to %s (%d kb windows, "
                         "%d um tiles, %s)\n"
                         % (len(windows), args.window_out,
                            args.window_size // 1000, args.window_tile_um,
                            ",".join(sorted(win_chroms))))

    if snp_counts is not None:
        # Every SNP that was ever informative, with the two alleles it was
        # declared with. The declared alleles are written out verbatim so a
        # suspected flip can be read straight off this file against the
        # reference base, without re-deriving the bed's column layout.
        n_snp = n_obs_tot = 0
        with gzip.open(args.snp_out, "wt") as sf:
            sf.write("chrom\tpos\tref_allele\talt_allele\tobs_ref\tobs_alt"
                     "\tmol_ref\tmol_alt\tgene\n")
            for c in sorted(snp_counts):
                smap = snps[c][1]
                # NOT `pos`: that name holds the barcode -> coordinate dict for
                # the whole of main(), and rebinding it here silently breaks
                # every use of it below.
                for sp in sorted(snp_counts[c]):
                    o_r, o_a, m_r, m_a = snp_counts[c][sp]
                    al = smap.get(sp)
                    g = gene_iv.name_at(c, sp) if gene_iv is not None else None
                    sf.write("%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n"
                             % (c, sp,
                                al[0] if al else ".", al[1] if al else ".",
                                o_r, o_a, m_r, m_a, g if g else ""))
                    n_snp += 1
                    n_obs_tot += o_r + o_a
        sys.stderr.write("Wrote %d informative SNPs (%d base observations) on "
                         "%s to %s\n"
                         % (n_snp, n_obs_tot, ",".join(sorted(snp_counts)),
                            args.snp_out))

    if gene_bins is not None:
        # Only on-tissue pixels, unlike the locus table: this one is a spatial
        # matrix and an unresolved bin has no pixel to occupy.
        with gzip.open(args.gene_bin_out, "wt") as gf:
            gf.write("chrom\tgene\trow\tcol\tref\talt\n")
            for (c, g, r_, c_), (r, a) in sorted(gene_bins.items()):
                gf.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (c, g, r_, c_, r, a))
        n_genes = len(set(k[1] for k in gene_bins))
        n_px = len(set((k[2], k[3]) for k in gene_bins))
        sys.stderr.write("Wrote %d gene x pixel rows (%d genes, %d pixels of "
                         "%dum) to %s\n"
                         % (len(gene_bins), n_genes, n_px, args.gene_bin_um,
                            args.gene_bin_out))

    if loci is not None:
        # NOTE: unlike the main table, this one counts every barcode the BAM
        # carried, including bins that did not resolve to tissue. That is
        # deliberate - its job is to check a control locus's depth and direction,
        # and the off-tissue bins carry the same alleles - but it means the totals
        # here run a few percent above the corresponding column in the main table.
        with open(args.locus_out, "w") as lf:
            lf.write("subset\tlocus\tref\talt\tn\tref_frac\n")
            for (sn, nm), (r, a) in sorted(loci.items(),
                                           key=lambda kv: (kv[0][0], -sum(kv[1]))):
                lf.write("%s\t%s\t%d\t%d\t%d\t%.4f\n"
                         % (sn, nm, r, a, r + a,
                            r / (r + a) if r + a else float("nan")))
        sys.stderr.write("Wrote %d subset x locus rows to %s\n"
                         % (len(loci), args.locus_out))

    n_bins = len(pos)
    match_rate = (100.0 * stats["snp_matching_an_allele"] / stats["snp_observations"]
                  if stats["snp_observations"] else 0.0)

    unit_lab = "UMIs" if args.count_unit == "molecule" else "reads"
    # Counted units per independent molecule - the number every SE downstream
    # has to be deflated by. In molecule mode the unit already IS the molecule,
    # so it is 1 by construction. In read mode it is measured, and it is NOT
    # reads-per-non-duplicate-read: dropping duplicates still leaves ~2.5 reads
    # per molecule here, because the duplicate flag is per alignment position
    # and a 3' molecule spans several. Deflating by the wrong one of those two
    # is what would let a read-level sweep recommend a tile size that cannot
    # measure what it claims to.
    n_mol = len(mol_seen) if mol_seen is not None else 0
    if args.count_unit == "molecule":
        dup_factor = 1.0
    elif n_mol:
        dup_factor = float(stats["informative_reads_primary"]) / n_mol
    else:
        dup_factor = float("nan")

    sys.stderr.write("\n--- summary ---\n")
    sys.stderr.write("counting unit              %s%s\n"
                     % (args.count_unit,
                        "" if args.drop_duplicates else ", PCR DUPLICATES KEPT"))
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
    sys.stderr.write("%-26s %d\n" % (unit_lab + " dropped as ties",
                                      stats["unit_ties"]))
    sys.stderr.write("informative reads          %d, of which %d survive the\n"
                     "  duplicate flag\n"
                     % (stats["informative_reads"],
                        stats["informative_reads_nondup"]))
    if args.count_unit == "read":
        sys.stderr.write("distinct (bin, UB) molecules %d -> duplication factor %.3f\n"
                         % (n_mol, dup_factor))
        sys.stderr.write(
            "  *** These counts are NOT %d independent observations. Divide n\n"
            "  *** by %.3f before quoting any binomial SE; ase_tile_sweep.R does\n"
            "  *** that from the duplication_factor in the provenance sidecar,\n"
            "  *** and reports the result as n_eff.\n"
            % (stats["informative_reads_primary"], dup_factor))
    sys.stderr.write("bins with informative %-5s %d / %d\n"
                     % (unit_lab, n_written, n_bins))
    sys.stderr.write("chrX informative %-10s %d   (lambda = %.4f per 2um bin)\n"
                     % (unit_lab, x_umi, x_umi / n_bins if n_bins else 0.0))
    sys.stderr.write("autosomal informative %-5s %d   (lambda = %.4f per 2um bin)\n"
                     % (unit_lab, a_umi, a_umi / n_bins if n_bins else 0.0))
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

    # Both from the write loop, so both count only bins that resolved to tissue.
    # Summing the numerator over every barcode in sub_counts while taking the
    # denominator from the write loop inflated every subset ref fraction by the
    # off-tissue share - 1.9% in 9w, 5.3% in 78w - which is enough to make the
    # subsets fail to add back up to the chrX column and to misreport the Xic.
    # The written table was never affected; only this line was.
    for i, s in enumerate(subsets):
        r, n = sub_ref[i], sub_umi[i]
        sys.stderr.write("subset %-12s %d informative %s, ref fraction %s\n"
                         % (s.name, n, unit_lab,
                            "%.4f" % (r / n) if n else "n/a"))

    # ------------------------------------------------------------- provenance
    #
    # Which SNP bed produced a result is not recoverable after the fact: three
    # "no-Xist" builds are named in this repo and only one is live, and none of
    # the output files carries a fingerprint. Written unconditionally, next to
    # the counts, so an archived result can always be attributed.
    prov = args.provenance_out or (args.out + ".provenance.tsv")
    with open(prov, "w") as pf:
        pf.write("key\tvalue\n")
        pf.write("argv\t%s\n" % " ".join(sys.argv))
        pf.write("bam\t%s\n" % args.bam)
        pf.write("spatial_dir\t%s\n" % args.spatial_dir)
        pf.write("snp_bed\t%s\n" % args.snps)
        pf.write("snp_bed_md5\t%s\n" % md5_of(args.snps))
        pf.write("snp_offset\t%d\n" % args.snp_offset)
        pf.write("min_mapq\t%d\n" % args.min_mapq)
        pf.write("min_baseq\t%d\n" % args.min_baseq)
        pf.write("drop_duplicates\t%s\n" % args.drop_duplicates)
        pf.write("count_unit\t%s\n" % args.count_unit)
        # Counted units per distinct molecule: 1.0 in molecule mode, measured in
        # read mode. ase_tile_sweep.R reads this key to turn a read-level n into
        # an effective one, so renaming it breaks the SEs downstream.
        pf.write("duplication_factor\t%.4f\n" % dup_factor)
        pf.write("informative_reads\t%d\n" % stats["informative_reads"])
        pf.write("informative_reads_nondup\t%d\n"
                 % stats["informative_reads_nondup"])
        pf.write("informative_molecules\t%d\n" % n_mol)
        # Task 9's rule: every bed that shaped a result is named and hashed
        # next to it, so a gene table can be attributed to its annotation.
        if args.gene_bin_out:
            pf.write("gene_bin_out\t%s\n" % args.gene_bin_out)
            pf.write("gene_bin_um\t%d\n" % args.gene_bin_um)
            pf.write("gene_assignment\t%s\n"
                     % ("gene_bed" if args.gene_bed else args.gene_tag + "_tag"))
            if args.gene_bed:
                pf.write("gene_bed\t%s\n" % args.gene_bed)
                pf.write("gene_bed_md5\t%s\n" % md5_of(args.gene_bed))
        pf.write("x_chrom\t%s\n" % args.x_chrom)
        pf.write("autosome_chroms\t%s\n" % ",".join(autosomes))
        for s in subsets:
            pf.write("subset_%s_bed\t%s\n" % (s.name, s.path))
            pf.write("subset_%s_bed_md5\t%s\n" % (s.name, md5_of(s.path)))
            pf.write("subset_%s_complement\t%s\n" % (s.name, s.complement))
            pf.write("subset_%s_intervals\t%d\n" % (s.name, s.iv.n))
        # How many distinct SNPs and genes actually contribute. Without this the
        # 12% pair-correlation floor cannot be attributed: one badly behaved
        # locus carrying a tenth of the chrX signal and a thousand loci each
        # carrying a thousandth look identical in the aggregate.
        for lab in sorted(seen):
            pf.write("informative_snps_%s\t%d\n" % (lab, len(seen[lab]["snps"])))
            pf.write("informative_genes_%s\t%d\n" % (lab, len(seen[lab]["genes"])))
        pf.write("informative_umis_chrX\t%d\n" % x_umi)
        pf.write("informative_umis_autosome\t%d\n" % a_umi)
        for i, s in enumerate(subsets):
            pf.write("informative_umis_%s\t%d\n" % (s.name, sub_umi[i]))
    sys.stderr.write("\nProvenance written to %s\n" % prov)
    for lab in sorted(seen):
        sys.stderr.write("  %-12s %d distinct informative SNPs, %d distinct %s "
                         "values\n" % (lab, len(seen[lab]["snps"]),
                                       len(seen[lab]["genes"]), args.gene_tag))
    sys.stderr.write("  (%s counts EXONIC assignment only, so it is a floor on "
                     "the gene count,\n   not the number of genes contributing "
                     "- most reads here are intronic.)\n" % args.gene_tag)

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
