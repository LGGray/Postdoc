#!/bin/bash
# Build ONE IGV session holding the four bulk-RNA conditions and the H3K27ac track,
# each as a B6/CAST overlay row, with every bigWig copied into a single self-contained
# folder so the session can be rsync'd to a laptop as a unit.
#
# Row order is the figure order: adult, Sham, TAC, aged, then H3K27ac.
#
#   bash slurm/build_figure_session.sh [LOCUS] [RNA_MAX] [CHIP_MAX] [ALSUB]
#
#     $1 LOCUS     region the session opens on, and the region the y-axes are scaled to
#                  (default chr5:148922383-148998282, the 5930430L01Rik / Hmgb1 locus)
#     $2 RNA_MAX   y-axis ceiling shared by ALL FOUR RNA rows   (default: auto)
#     $3 CHIP_MAX  y-axis ceiling for the H3K27ac row           (default: auto)
#     $4 ALSUB     which Allelome.LINK run to take the bed/bedpe from:
#                  AL | AL_65 | AL_chrX                         (default AL)
#
# Each condition contributes three rows, matching the layout of the Hmgb1 figure:
#   <label>_links.bedpe   Allelome.LINK arcs, red = repressing, green = enhancing
#   <label>_loci.bed      Allelome.PRO loci, coloured by allelic ratio, names carrying
#                         reads and ratio (Hmgb1_reads:143/36_score:9.847_ratio:0.77)
#   <label> coverage      the B6/CAST overlay
# H3K27ac has no links, so it contributes its peak bed and the overlay only.
#
# "auto" reads the actual maximum of each bigWig across LOCUS and rounds up to a round
# number. The four RNA rows always share one ceiling so their heights mean the same
# thing; H3K27ac gets its own, because ChIP coverage and RNA coverage are different
# quantities and forcing them onto one axis would flatten whichever is smaller.
#
# Cheap enough to run on a login node - it only reads bigWig headers and summaries -
# but it is equally happy inside the interactive allocation.

set -eo pipefail          # no -u: conda's activate.d hooks read unset variables

LOCUS=${1:-chr5:148922383-148998282}
RNA_MAX=${2:-auto}
CHIP_MAX=${3:-auto}
ALSUB=${4:-AL}

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

ROOT=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
RNADIR=$ROOT/heart_bulk_allelic_igv
CHIPDIR=$ROOT/heart_H3K27ac/igv
FIG=$ROOT/figure_igv

# label | bigWig dir | data dir (Allelome output lives here) | display name | group
ROWS=(
  "He_9w|$RNADIR|$ROOT/adult_aged_bodymap/He_9w|adult 9w   [BxC, random XCI]|rna"
  "He_Sham|$RNADIR|$ROOT/F1_TAC_Sarah/He_TACSham28d_RNA_XistBxC_-+_XX|Sham 28d  [XistBxC -+, CAST = Xi]|rna"
  "He_TAC|$RNADIR|$ROOT/F1_TAC_Sarah/He_TAC28d_RNA_XistBxC_-+_XX|TAC 28d   [XistBxC -+, CAST = Xi]|rna"
  "He_78w|$RNADIR|$ROOT/adult_aged_bodymap/He_78w|aged 78w  [BxC, random XCI]|rna"
  "H3K27ac|$CHIPDIR|$ROOT/heart_H3K27ac|H3K27ac  adult heart|chip"
)

mkdir -p "$FIG"

######------ collect the bigWigs into one folder ------######
missing=0
take () {                 # take <source> <name in $FIG>; records misses instead of dying
  if [[ -f $1 ]]; then
    cp -f "$1" "$FIG/$2"  # always refresh: a rerun of the track scripts changes these
  else
    echo "MISSING: $1" >&2
    missing=1
  fi
}

for entry in "${ROWS[@]}"; do
  IFS='|' read -r lab bwdir data _ grp <<< "$entry"
  for allele in B6 CAST; do
    take "$bwdir/${lab}_${allele}.bw" "${lab}_${allele}.bw"
  done
  if [[ $grp == chip ]]; then
    # the replicate-averaged peak table written by chipseq.slurm
    take "$data/H3K27ac_adult_heart.bed" "${lab}_loci.bed"
  else
    take "$data/$ALSUB/$ALSUB.bed"   "${lab}_loci.bed"
    take "$data/$ALSUB/$ALSUB.bedpe" "${lab}_links.bedpe"
  fi
done
if (( missing )); then
  echo "-> rows above will render empty until those bigWigs exist" >&2
fi

######------ pick the y-axis ceilings ------######
# pyBigWig ships with deeptools, so it is already in the RNAseq env.
locus_max () {          # locus_max <bigwig> [<bigwig> ...] -> max value over LOCUS
  python3 - "$LOCUS" "$@" <<'PY'
import sys
try:
    import pyBigWig
except ImportError:
    print("nan"); sys.exit(0)
locus = sys.argv[1]
chrom, rng = locus.split(":")
start, end = (int(x.replace(",", "")) for x in rng.split("-"))
best = 0.0
for path in sys.argv[2:]:
    try:
        bw = pyBigWig.open(path)
        v = bw.stats(chrom, start, end, type="max")[0]
        bw.close()
        if v:
            best = max(best, float(v))
    except Exception as e:
        print(f"could not read {path}: {e}", file=sys.stderr)
print(f"{best:.6f}")
PY
}

# round up to the next 1/2/5 x 10^n, so the axis label stays readable
nice_ceiling () {
  awk -v v="$1" 'BEGIN{
    if (v != v || v <= 0) { print "nan"; exit }         # nan guard
    p = 1; while (v/p >= 10) p *= 10; while (v/p < 1) p /= 10;
    m = v/p;
    c = (m <= 1) ? 1 : (m <= 2) ? 2 : (m <= 5) ? 5 : 10;
    printf "%g", c*p
  }'
}

if [[ $RNA_MAX == auto ]]; then
  raw=$(locus_max "$FIG"/He_9w_*.bw "$FIG"/He_78w_*.bw "$FIG"/He_Sham_*.bw "$FIG"/He_TAC_*.bw)
  RNA_MAX=$(nice_ceiling "$raw")
  if [[ $RNA_MAX == nan ]]; then RNA_MAX=20; echo "auto-scale failed for RNA, using 20" >&2; fi
  echo "RNA  max over $LOCUS: $raw -> axis 0-$RNA_MAX"
fi
if [[ $CHIP_MAX == auto ]]; then
  raw=$(locus_max "$FIG"/H3K27ac_*.bw)
  CHIP_MAX=$(nice_ceiling "$raw")
  if [[ $CHIP_MAX == nan ]]; then CHIP_MAX=10; echo "auto-scale failed for ChIP, using 10" >&2; fi
  echo "ChIP max over $LOCUS: $raw -> axis 0-$CHIP_MAX"
fi

######------ write the session ------######
IGVLOCUS=${LOCUS//,/}
{
  cat <<XML
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" locus="${IGVLOCUS}" version="8">
  <Resources>
XML
  for entry in "${ROWS[@]}"; do
    IFS='|' read -r lab _ _ _ grp <<< "$entry"
    printf '    <Resource path="%s_B6.bw"/>\n    <Resource path="%s_CAST.bw"/>\n' "$lab" "$lab"
    printf '    <Resource path="%s_loci.bed"/>\n' "$lab"
    if [[ $grp != chip ]]; then printf '    <Resource path="%s_links.bedpe"/>\n' "$lab"; fi
  done
  echo '  </Resources>'
  echo '  <Panel name="DataPanel">'
  for entry in "${ROWS[@]}"; do
    IFS='|' read -r lab _ _ disp grp <<< "$entry"
    if [[ $grp == chip ]]; then max=$CHIP_MAX; else max=$RNA_MAX; fi
    # arcs on top, then the loci, then the coverage overlay - the figure's order
    if [[ $grp != chip ]]; then
      cat <<XML
    <Track clazz="org.broad.igv.track.FeatureTrack" id="${lab}_links.bedpe" attributeKey="${lab}_links"
           name="${disp}  links" height="60" visible="true" displayMode="COLLAPSED"/>
XML
    fi
    cat <<XML
    <Track clazz="org.broad.igv.track.FeatureTrack" id="${lab}_loci.bed" attributeKey="${lab}_loci"
           name="${disp}  loci" height="25" visible="true" displayMode="COLLAPSED"/>
XML
    cat <<XML
    <Track clazz="org.broad.igv.track.MergedTracks" id="${lab}_allelic" attributeKey="${lab}"
           name="${disp}" alpha="0.5" autoScale="false"
           renderer="BAR_CHART" height="80" visible="true">
      <DataRange minimum="0.0" baseline="0.0" maximum="${max}" type="LINEAR"/>
      <Track clazz="org.broad.igv.track.DataSourceTrack" id="${lab}_B6.bw" attributeKey="${lab}_B6.bw"
             name="${lab} B6" color="0,90,181" autoScale="false"
             renderer="BAR_CHART" windowFunction="mean" visible="true">
        <DataRange minimum="0.0" baseline="0.0" maximum="${max}" type="LINEAR"/>
      </Track>
      <Track clazz="org.broad.igv.track.DataSourceTrack" id="${lab}_CAST.bw" attributeKey="${lab}_CAST.bw"
             name="${lab} CAST" color="220,50,32" autoScale="false"
             renderer="BAR_CHART" windowFunction="mean" visible="true">
        <DataRange minimum="0.0" baseline="0.0" maximum="${max}" type="LINEAR"/>
      </Track>
    </Track>
XML
  done
  cat <<XML
  </Panel>
  <Panel name="FeaturePanel">
    <Track clazz="org.broad.igv.track.SequenceTrack" id="Reference sequence" name="Reference sequence"/>
  </Panel>
</Session>
XML
} > "$FIG/.figure_session.xml.tmp"
mv "$FIG/.figure_session.xml.tmp" "$FIG/figure_session.xml"

echo
echo "session: $FIG/figure_session.xml"
echo "RNA rows 0-$RNA_MAX, H3K27ac row 0-$CHIP_MAX, opening on $IGVLOCUS"
ls -lh "$FIG"
echo
echo "copy to the laptop with:"
echo "  rsync -av <cluster>:$FIG/ ~/igv/figure/"
