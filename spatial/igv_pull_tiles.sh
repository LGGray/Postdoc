#!/usr/bin/env bash
set -euo pipefail
trap 'st=$?; echo "ERROR: $0 stopped at line $LINENO (exit $st)" >&2' ERR

# Pull a handful of tile BAMs into one directory for IGV: N CAST-biased ("blue")
# tiles and N depth-matched strongly-Bl6 ("red") ones per sample, each in BOTH
# the duplicate-inclusive and the deduplicated form.
#
# WHY BOTH FORMS, and this is the whole point of the exercise:
#
#   IGV HIDES DUPLICATE-FLAGGED READS BY DEFAULT
#   (View > Preferences > Alignments > "Filter duplicate reads")
#
# which is exactly what `samtools mpileup` does inside Allelome.PRO2. So with
# IGV's default settings the _flt file and the _dup file look the SAME, and both
# show what Allelome.PRO2 scored. Untick that preference, or right-click the
# track and turn the duplicate filter off, and the _dup file shows what is
# actually in the library. The blue tiles are blue only in the filtered view:
# 84.5% of informative chrX reads in this data are duplicate-flagged and ~95% of
# those are Bl6, so removing them strips the Bl6 evidence and the tile flips.
# Measured on the 421 CAST-skewed tiles: every one is Bl6-dominant at molecule
# level (78w median 0.818 against Allelome.PRO2's 0.229).
#
# So the demonstration for the boss is: load a blue tile, look at it with the
# duplicate filter ON (blue, as Allelome.PRO2 sees it), then OFF (the pile-ups
# appear and it is not blue). A red tile is the control - it looks the same
# either way.
#
# USAGE (login node; this is a copy, not a computation)
#   conda activate RNAseq                      # only needed if a .bai is missing
#   ./spatial/igv_pull_tiles.sh                        # 3+3 per sample, both samples
#   ./spatial/igv_pull_tiles.sh ~/igv_tiles 3          # explicit outdir and N
#   SAMPLES=78w ./spatial/igv_pull_tiles.sh            # exactly 3 blue + 3 red
#   CSV=/path/to/other_map.csv ./spatial/igv_pull_tiles.sh
#
# Writes into OUTDIR: the BAMs and indexes, manifest.tsv (what each file is and
# why it was picked), and loci.txt (the chrX genes carrying each tile's reads,
# so IGV has somewhere to jump to - a tile BAM spans the whole genome).

BASE=${BASE:-/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2}
SPATIAL=$BASE/adult_aged_spatial

OUTDIR=${1:-$SPATIAL/ase/igv_tiles}
N=${2:-3}
SAMPLES=${SAMPLES:-9w,78w}
TILE_UM=${TILE_UM:-64}
# The floored map's CSV, because its `call` column is the one that says which
# tiles Allelome.PRO2 called CAST-skewed against this sample's own null band.
CSV=${CSV:-$SPATIAL/ase/tile_ratio_map_${TILE_UM}um_dup_min20.csv}
# "High red" = the top OCM bins. 0.95 is the boundary of the reddest bin.
RED_MIN=${RED_MIN:-0.95}

if [ ! -f "$CSV" ]; then
  echo "ERROR: no map CSV at $CSV" >&2
  echo "       run slurm/spatial_tile_floor_map.slurm first, or set CSV=" >&2
  exit 1
fi
if ! [[ "$N" =~ ^[0-9]+$ ]] || [ "$N" -lt 1 ]; then
  echo "ERROR: N='$N' must be a positive integer" >&2; exit 1
fi

# Column positions, resolved from the header rather than hardcoded: the floor
# script appends columns (lowratio) and a fixed $16 would silently read the
# wrong field the next time one is added.
# `|| true` so a missing column yields an empty string and the explicit check
# below reports which one, rather than `set -e` killing the assignment first.
col() { head -1 "$CSV" | tr ',' '\n' | grep -nxF "$1" | cut -d: -f1 || true; }
C_SAMPLE=$(col sample); C_TILE=$(col tile); C_A1=$(col x_a1); C_A2=$(col x_a2)
C_N=$(col x_n); C_RATIO=$(col x_ratio); C_Z=$(col z); C_CALL=$(col call)
for v in C_SAMPLE C_TILE C_A1 C_A2 C_N C_RATIO C_Z C_CALL; do
  [ -n "${!v}" ] || { echo "ERROR: $CSV has no column for $v" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
MAN=$OUTDIR/manifest.tsv
LOCI=$OUTDIR/loci.txt
printf 'sample\tclass\ttile\tchrX_reads\tA1_Bl6\tA2_CAST\tratio\tz\tfiles\n' > "$MAN"
: > "$LOCI"

echo "map CSV:   $CSV"
echo "out dir:   $OUTDIR"
echo "picking:   $N CAST-biased + $N depth-matched Bl6 tiles per sample"
echo

# Copy one tile BAM plus its index, from the first root that has it. bam_dup is
# scratch-only; bam_flt exists in both places and permanent storage is
# preferred, since $SCRATCH is purged on an access-time policy.
copy_tile() {
  local smp=$1 tile=$2 sub=$3 suffix=$4 dest=$5 root found=""
  for root in "$SPATIAL/ase/$smp/$sub" \
              "${SCRATCH:-/nonexistent}/spatial_tiles_${smp}_${TILE_UM}um/$sub"; do
    [ -f "$root/$tile.bam" ] && { found=$root; break; }
  done
  if [ -z "$found" ]; then
    echo "    MISSING $sub/$tile.bam (looked in ase/$smp and \$SCRATCH)" >&2
    return 1
  fi
  cp "$found/$tile.bam" "$dest/${smp}_${tile}${suffix}.bam"
  if [ -f "$found/$tile.bam.bai" ]; then
    cp "$found/$tile.bam.bai" "$dest/${smp}_${tile}${suffix}.bam.bai"
  else
    # IGV will not load an unindexed BAM. Index the copy rather than the
    # source, which may be read-only or on scratch.
    echo "    no .bai for $sub/$tile - indexing the copy"
    samtools index "$dest/${smp}_${tile}${suffix}.bam"
  fi
  echo "${smp}_${tile}${suffix}.bam"
}

# The chrX genes actually carrying the tile's reads, so IGV has a locus list.
# Best effort: the gene-level table is the pysam one, and the run is fine
# without it - IGV just opens on whatever locus it opened on last.
GENE_TSV_FOUND=0
top_genes() {
  local smp=$1 tile=$2 gz=$SPATIAL/ase_pysam_dup_${TILE_UM}um/$smp/tile_gene_counts.tsv.gz
  [ -f "$gz" ] || return 0
  GENE_TSV_FOUND=1
  gzip -dc "$gz" | awk -F'\t' -v t="$tile" -v s="$smp" '
    NR>1 && $1==t && $4=="chrX" && $11+0 > 0 {
      printf "%s\t%s\t%s\t%s\tA1=%s A2=%s reads=%s\n", $11, s, t, $5, $9, $10, $11 }' \
    | sort -k1,1nr | awk 'NR <= 5' | cut -f2-
}

# Defined once, not per sample: $smp and $OUTDIR resolve when it is called.
emit() {                        # class, then the tab-separated tile record
  local cls=$1 n=$2 tile=$3 a1=$4 a2=$5 ratio=$6 z=$7 files=""
  # 15 significant figures of a ratio is noise the reader has to step over.
  ratio=$(printf '%.3f' "$ratio"); z=$(printf '%.2f' "$z")
  printf '  %-11s %-24s %5s reads  A1=%-6s A2=%-6s ratio=%s\n' \
         "$cls" "$tile" "$n" "$a1" "$a2" "$ratio"
  for pair in "bam_dup:_dup" "bam_flt:_flt"; do
    local sub=${pair%%:*} sfx=${pair##*:} f
    if f=$(copy_tile "$smp" "$tile" "$sub" "$sfx" "$OUTDIR"); then
      files="${files:+$files,}$f"
    fi
  done
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
         "$smp" "$cls" "$tile" "$n" "$a1" "$a2" "$ratio" "$z" "$files" >> "$MAN"
  top_genes "$smp" "$tile" >> "$LOCI"
}


for smp in ${SAMPLES//,/ }; do
  echo "[$smp]"

  blue=$(awk -F, -v s="$smp" -v cs=$C_SAMPLE -v cc=$C_CALL -v cn=$C_N -v ct=$C_TILE \
             -v c1=$C_A1 -v c2=$C_A2 -v cr=$C_RATIO -v cz=$C_Z \
    'NR>1 && $cs==s && $cc=="CAST-skewed" {print $cn"\t"$ct"\t"$c1"\t"$c2"\t"$cr"\t"$cz}' \
    "$CSV" | sort -k1,1nr | awk -v n="$N" 'NR <= n')
  if [ -z "$blue" ]; then
    echo "  no CAST-skewed tiles in $CSV for $smp - skipping"; continue
  fi

  # Depth-match the red tiles to the blue ones. Without this the comparison is
  # confounded: the CAST-skewed tiles are the DEEPEST on the slide (78w median
  # 136 chrX reads against 53 for the Bl6-skewed ones), so the N deepest red
  # tiles would not be a like-for-like control.
  med=$(echo "$blue" | cut -f1 | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
  echo "  blue: $N deepest CAST-skewed tiles, median depth $med"
  red=$(awk -F, -v s="$smp" -v cs=$C_SAMPLE -v cr=$C_RATIO -v cn=$C_N -v ct=$C_TILE \
            -v c1=$C_A1 -v c2=$C_A2 -v cz=$C_Z -v rmin="$RED_MIN" -v m="$med" \
    'NR>1 && $cs==s && $cr+0 >= rmin {d=$cn-m; if (d<0) d=-d;
     print d"\t"$cn"\t"$ct"\t"$c1"\t"$c2"\t"$cr"\t"$cz}' \
    "$CSV" | sort -k1,1n | awk -v n="$N" 'NR <= n' | cut -f2-)
  if [ -z "$red" ]; then
    echo "  WARNING: no tile reaches ratio >= $RED_MIN in $smp - lower RED_MIN" >&2
  fi
  echo "  red:  $N tiles at ratio >= $RED_MIN, closest in depth to $med"

  while IFS=$'\t' read -r n tile a1 a2 ratio z; do
    [ -n "${tile:-}" ] && emit CAST-biased "$n" "$tile" "$a1" "$a2" "$ratio" "$z"
  done <<< "$blue"
  while IFS=$'\t' read -r n tile a1 a2 ratio z; do
    [ -n "${tile:-}" ] && emit Bl6-biased "$n" "$tile" "$a1" "$a2" "$ratio" "$z"
  done <<< "$red"
  echo
done

echo "Wrote $(ls "$OUTDIR"/*.bam 2>/dev/null | wc -l | tr -d ' ') BAMs to $OUTDIR"
echo "      $MAN"
if [ "$GENE_TSV_FOUND" -eq 1 ]; then
  echo "      $LOCI   (sample, tile, gene, counts - search the gene name in IGV)"
else
  rm -f "$LOCI"
  echo "      no ase_pysam_dup_${TILE_UM}um gene table found, so no locus list"
fi
cat <<'NOTE'

IN IGV, and this is the part that matters:

  1. Load a _dup BAM and its _flt twin for the same tile.
  2. With IGV's default settings they look IDENTICAL - IGV filters
     duplicate-flagged reads, exactly as Allelome.PRO2's mpileup does. This is
     the view in which the tile is blue.
  3. View > Preferences > Alignments > untick "Filter duplicate reads"
     (or right-click the track > "Duplicates" toggle, depending on version).
     The _dup track fills in with pile-ups; the _flt track does not change,
     because those reads were removed from the file itself.
  4. Colour the track by "read strand" off and use "Shade base by quality" /
     "Colour alignments by > tag" if you want the allele visible per read;
     the SNP columns are what carry it.

A CAST-biased tile is blue in step 2 and not in step 3. A Bl6-biased tile
looks the same in both, which is why the red tiles are here as controls.

Start from loci.txt rather than browsing: a tile BAM spans the whole genome and
the chrX reads are a thin scatter across it. The genes listed there are where
each tile's reads actually are, and the ones supplying the A2 (CAST) reads are
worth checking against ase/<sample>/artifact_snps_<sample>.tsv - a locus that
is near-pure CAST across the whole section is a SNP problem, not a tile.
NOTE
