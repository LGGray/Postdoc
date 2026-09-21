#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J chip_AS_track
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8    # serial_std caps at 16, but 8 is the real need: haplotag is
                             # single-threaded per contig and 8 workers already cover the
                             # makespan, which is bound by chr1 alone
#SBATCH --mem=24G
#SBATCH --time=8:00:00       # serial_std hard cap is 24h; the run is ~2-4h, and a smaller
                             # ask backfills far sooner

# Allele-specific ChIP-seq tracks for IGV.
#
#   merged.clean.sorted.bam
#     -> whatshap haplotag  (B6/CAST SNPs as a pseudo-phased VCF, one chromosome per core)
#     -> whatshap split     -> <MARK>_B6.bam   (H1 = REF = B6   = maternal)
#                           -> <MARK>_CAST.bam (H2 = ALT = CAST = paternal)
#     -> bamCoverage        -> <MARK>_B6.bw    (blue)
#                           -> <MARK>_CAST.bw  (red), same scale factor
#
# The session XML overlays the two into one IGV row (a MergedTracks element, which IGV
# draws at 50% alpha): blue where B6 dominates, red where CAST does, purple where they
# match. Both bigWigs are therefore POSITIVE - do not negate CAST, or the overlay puts
# one allele off the bottom of the panel.
#
# Both bigWigs use ONE shared scale factor, 1e6 / (N_B6 + N_CAST). Normalising each
# track to its own library size would divide the allelic imbalance away, which is the
# only thing the track is meant to show. With a shared denominator the two bar heights
# are directly comparable and their sum is the assignable coverage in CPM.
#
#
# ---- TWO WAYS TO RUN THIS -------------------------------------------------------
#
# (a) batch, --export=NONE so arguments are positional:
#
#       sbatch slurm/chipseq_allelic_tracks.slurm H3K27ac
#
#     $1 MARK     H3K27ac | H3K4me3            (default H3K27ac -> heart_H3K27ac/)
#     $2 FRAGLEN  read extension in bp         (default 188; single-end ChIP, so bamCoverage
#                                               needs an explicit fragment length. 188 is the
#                                               mean of the MACS2 "alternative fragment length"
#                                               for the two H3K27ac replicates, 183 and 193.
#                                               H3K4me3: check its own *_peaks.xls.)
#     $3 BINSIZE  bigWig bin size in bp        (default 25)
#     $4 MAXV     y-axis limit in the IGV session XML (default 10)
#
# (b) interactively, block by block. Grab an allocation first:
#
#       salloc --clusters=inter --partition=cm4_inter --nodes=1 \
#              --cpus-per-task=8 --mem=24G --time=8:00:00
#       srun --pty bash -l                  # shell ON the allocated node
#
#     then select and run BLOCK 0 (it defines every variable the later blocks use),
#     then BLOCK 1..5 in order. Blocks are re-runnable: BLOCK 1 skips an existing VCF
#     and BLOCK 2 skips contigs it has already tagged, so a block that dies part-way
#     can simply be run again. To override a default interactively, edit the
#     assignment in BLOCK 0 before running it - there is no $1 in a pasted shell.
#
#     BLOCK 2 is the long one (~45-60 min) and BLOCK 3 is single-threaded. If the
#     session might drop, run those inside tmux.
# ---------------------------------------------------------------------------------

# Strict mode under sbatch only. Pasted into an interactive shell, `set -e` would close
# the shell on the first failed command, and `set -u` on any unset variable.
case $- in *i*) : ;; *) set -euo pipefail ;; esac


# ==================================================================================
# BLOCK 0  -  setup: paths, tools, defaults. Run this first, every session.
# ==================================================================================
date

MARK=${1:-H3K27ac}
FRAGLEN=${2:-188}
BINSIZE=${3:-25}
MAXV=${4:-10}

module load slurm_setup 2>/dev/null || true      # already loaded on an interactive node
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

ROOT=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
WORK=$ROOT/heart_${MARK}
BAM=$WORK/merged.clean.sorted.bam                       # replicate-merged, -F 1804 -q 30, blacklist-filtered
REF=$ROOT/bwa-index/GRCm39.primary_assembly.genome.fa
SNPBED=$ROOT/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed
VCF=$ROOT/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.vcf.gz
OUT=$WORK/igv
TMP=${SCRATCH:-$WORK}/whatshap_${MARK}                  # no job id in the name, so a rerun
                                                        # resumes from the same scratch dir

NPROC=${SLURM_CPUS_PER_TASK:-8}
JOBS=$(( NPROC < 8 ? NPROC : 8 ))                       # haplotag processes in flight; more
                                                        # than 8 buys little with 20 contigs
THRD_PER_JOB=$(( NPROC / JOBS )); THRD_PER_JOB=$(( THRD_PER_JOB > 0 ? THRD_PER_JOB : 1 ))

mkdir -p "$OUT" "$TMP"

# whatshap, samtools, htslib and deeptools all live in RNAseq, so no env switching.
conda activate RNAseq
for tool in whatshap samtools bgzip tabix bamCoverage; do
  command -v "$tool" >/dev/null || echo "MISSING TOOL: $tool is not in the RNAseq env" >&2
done
for f in "$BAM" "$REF" "$SNPBED"; do
  [[ -f $f ]] || echo "MISSING INPUT: $f" >&2
done

echo "mark    $MARK  ($BAM)"
echo "cores   $NPROC -> $JOBS haplotag workers x $THRD_PER_JOB output threads"
echo "scratch $TMP"
echo "output  $OUT"
echo "whatshap $(whatshap --version 2>/dev/null)"


# ==================================================================================
# BLOCK 1  -  SNP bed -> pseudo-phased VCF. Built once, then reused by every mark.
# ==================================================================================
# SNPfile_C57BL6_NJ_CAST_EiJ.bed is chr / pos / pos+1 / <REF><ALT> / 0 / .
# Column 2 is ALREADY the 1-based position, not a 0-based BED start: checked against
# GRCm39.primary_assembly.genome.fa, where 158617 of the 158624 chr1 SNPs in the first
# 20 Mb have their REF base at col2 and only 2 at col2+1. Do not "fix" this to $3.
# Column 4 is B6 base then CAST base, matching Allelome.PRO2's A1 = B6 convention, so
# GT 0|1 makes haplotype 1 = B6 and haplotype 2 = CAST.
if [[ -f ${VCF}.tbi ]]; then
  echo "VCF already built: $VCF"
else
  echo "building $VCF"
  {
    echo '##fileformat=VCFv4.2'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    # single pass, unlike the column-by-column paste in AS_BW_whatshap.sh: rows whose
    # 4th field is not two bases are dropped whole, so chrom/pos can never slip out of
    # register against REF/ALT
    awk 'BEGIN{OFS="\t"} $4 ~ /^[ACGT][ACGT]$/ {
           print $1, $2, ".", substr($4,1,1), substr($4,2,1), "50.0", ".", ".", "GT", "0|1"
         }' "$SNPBED"
  } | bgzip -c > "$VCF"
  tabix -p vcf "$VCF"
fi

[[ -f ${REF}.fai ]] || samtools faidx "$REF"
[[ -f ${BAM}.bai ]] || samtools index -@ "$NPROC" "$BAM"
echo "VCF contigs: $(tabix -l "$VCF" | tr '\n' ' ')"


# ==================================================================================
# BLOCK 2  -  haplotag, one chromosome per worker.  THE LONG ONE, ~45-60 min.
# ==================================================================================
# --regions makes haplotag fetch only that contig, so each worker reads its own slice
# of the BAM and writes its own read->haplotype list. The tagged BAMs are throwaway:
# only the lists are used downstream. Each list is written to .part and renamed on
# success, so an interrupted contig is redone rather than silently half-read.
tag_chrom () {
  chrom="$1"
  list="$TMP/${chrom}.haplotag.tsv"
  whatshap haplotag \
    --reference "$REF" \
    --regions "$chrom" \
    --ignore-read-groups --ignore-linked-read --skip-missing-contigs \
    --output-threads "$THRD_PER_JOB" \
    --output-haplotag-list "${list}.part" \
    -o "$TMP/${chrom}.tagged.bam" \
    "$VCF" "$BAM" || { echo "  FAILED $chrom" >&2; return 1; }
  rm -f "$TMP/${chrom}.tagged.bam"
  mv "${list}.part" "$list"
  echo "  done $chrom"
}
export -f tag_chrom
export REF VCF BAM TMP THRD_PER_JOB

mapfile -t CHROMS < <(tabix -l "$VCF")
echo "haplotagging ${#CHROMS[@]} contigs, $JOBS at a time"

# skip messages go to stderr: anything on stdout here is fed to xargs as a contig name
{ for chrom in "${CHROMS[@]}"; do
    if [[ -s $TMP/${chrom}.haplotag.tsv ]]; then
      echo "  skip $chrom (list already exists)" >&2
    else
      echo "$chrom"
    fi
  done
} | xargs -n1 -P "$JOBS" -I{} bash -c 'tag_chrom "$@"' _ {}

# one list for the whole genome: keep the first header, drop the rest
{ head -n1 "$TMP/${CHROMS[0]}.haplotag.tsv"
  cat "$TMP"/*.haplotag.tsv | grep -v '^#'
} > "$TMP/haplotag_list.tsv"
echo "tagged reads: $(grep -vc '^#' "$TMP/haplotag_list.tsv") in list, \
H1 $(awk -F'\t' '$2=="H1"' "$TMP/haplotag_list.tsv" | wc -l), \
H2 $(awk -F'\t' '$2=="H2"' "$TMP/haplotag_list.tsv" | wc -l)"


# ==================================================================================
# BLOCK 3  -  split into the two alleles, sort, index, count.  single-threaded, ~1h.
# ==================================================================================
# H1 = REF = B6 = maternal (the active X in this cross), H2 = ALT = CAST = paternal (Xi).
# Reads that never overlap a SNP are unassignable and are dropped here by design.
whatshap split \
  --output-h1 "$TMP/${MARK}_B6.bam" \
  --output-h2 "$TMP/${MARK}_CAST.bam" \
  "$BAM" "$TMP/haplotag_list.tsv"

for allele in B6 CAST; do
  samtools sort -@ "$NPROC" -m 1G -o "$OUT/${MARK}_${allele}.bam" "$TMP/${MARK}_${allele}.bam"
  samtools index -@ "$NPROC" "$OUT/${MARK}_${allele}.bam"
  rm -f "$TMP/${MARK}_${allele}.bam"
done

N_B6=$(samtools view -c -@ "$NPROC" "$OUT/${MARK}_B6.bam")
N_CAST=$(samtools view -c -@ "$NPROC" "$OUT/${MARK}_CAST.bam")
SF=$(awk -v a="$N_B6" -v b="$N_CAST" 'BEGIN{ printf "%.10f", 1000000/(a+b) }')
echo "assigned reads: B6 $N_B6, CAST $N_CAST"
awk -v a="$N_B6" -v b="$N_CAST" 'BEGIN{ printf "genome-wide B6 fraction: %.4f (escape on chrX is 1 - this, per peak)\n", a/(a+b) }'
echo "shared scale factor: $SF"


# ==================================================================================
# BLOCK 4  -  bigWigs.  ~20 min.  (needs $SF from BLOCK 3)
# ==================================================================================
bamCoverage -b "$OUT/${MARK}_B6.bam"   -o "$OUT/${MARK}_B6.bw" \
  --binSize "$BINSIZE" --extendReads "$FRAGLEN" --scaleFactor "$SF"  -p "$NPROC"
bamCoverage -b "$OUT/${MARK}_CAST.bam" -o "$OUT/${MARK}_CAST.bw" \
  --binSize "$BINSIZE" --extendReads "$FRAGLEN" --scaleFactor "$SF"  -p "$NPROC"


# ==================================================================================
# BLOCK 5  -  IGV session and summary.
# ==================================================================================
# Relative paths: copy the whole igv/ folder to the Mac and open the .xml in IGV.
# The two bigWigs are nested inside one MergedTracks element, which is what IGV writes
# when you pick "Overlay Tracks" in the GUI. Both tracks and the container carry the
# same DataRange: an overlay is only readable if the panes share one y-axis.
cat > "$OUT/${MARK}_allelic_session.xml" <<XML
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" locus="chrX:102503978-102526839" version="8">
  <Resources>
    <Resource path="${MARK}_B6.bw"/>
    <Resource path="${MARK}_CAST.bw"/>
  </Resources>
  <Panel name="DataPanel">
    <Track clazz="org.broad.igv.track.MergedTracks" id="${MARK}_allelic" attributeKey="${MARK} allelic"
           name="${MARK}  B6 blue / CAST red" alpha="0.5" autoScale="false"
           renderer="BAR_CHART" height="90" visible="true">
      <DataRange minimum="0.0" baseline="0.0" maximum="${MAXV}.0" type="LINEAR"/>
      <Track clazz="org.broad.igv.track.DataSourceTrack" id="${MARK}_B6.bw" attributeKey="${MARK}_B6.bw"
             name="B6 (maternal, Xa)" color="0,90,181" autoScale="false"
             renderer="BAR_CHART" windowFunction="mean" visible="true">
        <DataRange minimum="0.0" baseline="0.0" maximum="${MAXV}.0" type="LINEAR"/>
      </Track>
      <Track clazz="org.broad.igv.track.DataSourceTrack" id="${MARK}_CAST.bw" attributeKey="${MARK}_CAST.bw"
             name="CAST (paternal, Xi)" color="220,50,32" autoScale="false"
             renderer="BAR_CHART" windowFunction="mean" visible="true">
        <DataRange minimum="0.0" baseline="0.0" maximum="${MAXV}.0" type="LINEAR"/>
      </Track>
    </Track>
  </Panel>
  <Panel name="FeaturePanel">
    <Track clazz="org.broad.igv.track.SequenceTrack" id="Reference sequence" name="Reference sequence"/>
  </Panel>
</Session>
XML

echo "Wrote:"
ls -lh "$OUT"
echo "scratch left in place for reruns (delete when happy): $TMP"
date
