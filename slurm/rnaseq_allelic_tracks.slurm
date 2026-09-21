#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J rna_AS_track
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8    # serial_std caps at 16; haplotag is single-threaded per contig
#SBATCH --mem=24G
#SBATCH --time=8:00:00
#SBATCH --array=1-4%2        # one task per condition, two at a time
#
# Allele-specific bulk RNA tracks for IGV - the RNA counterpart of
# chipseq_allelic_tracks.slurm, same whatshap route and the same overlay output.
#
#   pooled.sorted.bam
#     -> whatshap haplotag  (B6/CAST SNPs as a pseudo-phased VCF, one chromosome per core)
#     -> whatshap split     -> <COND>_B6.bam   (H1 = REF = B6   = maternal)
#                           -> <COND>_CAST.bam (H2 = ALT = CAST = paternal)
#     -> bamCoverage        -> <COND>_B6.bw / <COND>_CAST.bw, both positive, one shared
#                              scale factor 1e6/(N_B6+N_CAST) per condition
#
# All four conditions write into ONE directory so a single session XML can stack them,
# and BLOCK 5 emits that session with one overlay row per condition.
#
#
# ---- READ THIS BEFORE INTERPRETING THE TRACKS -----------------------------------
#
# The four conditions are NOT one mouse model. From the sample names in
# adult_aged_bodymap/{adult,aged}_metadata.csv and F1_TAC_Sarah/metadata.csv:
#
#   He_9w    He_9w_RNA_BxC_XX_{1,2,3}                 BxC        wildtype F1
#   He_78w   He_78w_RNA_BxC_XX_{1,2,3}                BxC        wildtype F1
#   He_Sham  He_TACSham28d_RNA_XistBxC_-+_XX_{1,2,3}  XistBxC -+ Xist deleted on B6
#   He_TAC   He_TAC28d_RNA_XistBxC_-+_XX_{1,2,3}      XistBxC -+ Xist deleted on B6
#
# So the red (CAST) bar means two different things:
#
#   Sham/TAC   XCI is skewed - CAST is the inactive X in every cell, so CAST signal on
#              chrX is escape, exactly as in the H3K27ac track.
#   9w/78w     XCI is random - each cell silenced a different X, and a whole heart is a
#              mosaic. chrX sits near 50:50 whatever the escape level is, so a red bar
#              here is NOT escape and the two pairs must not be read off against each
#              other. Xist itself is the giveaway: skewed in Sham/TAC, biallelic-looking
#              in 9w/78w.
#
# The genotype is written into every track name so this cannot be lost in the figure.
# If the 9w/78w whole hearts are meant to be the Xist-deletion cross, they are the wrong
# directories - the XBxC -+ adult and aged hearts are the FACS sets in cardiac_RNAseq/
# and aged_cardiac_RNAseq/, and CONDITIONS below should be repointed.
# ---------------------------------------------------------------------------------
#
#
# ---- TWO WAYS TO RUN THIS -------------------------------------------------------
#
# (a) batch, --export=NONE so arguments are positional:
#
#       sbatch slurm/rnaseq_allelic_tracks.slurm
#
#     $1 BINSIZE  bigWig bin size in bp                   (default 10)
#     $2 MAXV     y-axis limit in the IGV session XML     (default 20)
#
# (b) interactively, block by block:
#
#       salloc --clusters=inter --partition=cm4_inter --nodes=1 \
#              --cpus-per-task=8 --mem=24G --time=8:00:00
#       srun --pty bash -l
#
#     Set IDX in BLOCK 0 to pick the condition (1..4), run BLOCK 0-4, then come back
#     and repeat with the next IDX. Run BLOCK 5 once at the end. Blocks are
#     re-runnable: BLOCK 1 skips an existing VCF, BLOCK 2 skips contigs already tagged.
#     Run BLOCK 2 and 3 inside tmux - together they are ~2h per condition.
# ---------------------------------------------------------------------------------

# Strict mode under sbatch only; see chipseq_allelic_tracks.slurm.
case $- in *i*) : ;; *) set -euo pipefail ;; esac


# ==================================================================================
# BLOCK 0  -  setup: condition table, paths, tools. Run this first, every session.
# ==================================================================================
date

BINSIZE=${1:-10}          # RNA is exon-sharp, so a finer bin than the ChIP track
MAXV=${2:-20}

module load slurm_setup 2>/dev/null || true
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

ROOT=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2

# label | directory | genotype shown in the track name
CONDITIONS=(
  "He_9w|$ROOT/adult_aged_bodymap/He_9w|BxC wildtype, random XCI"
  "He_78w|$ROOT/adult_aged_bodymap/He_78w|BxC wildtype, random XCI"
  "He_Sham|$ROOT/F1_TAC_Sarah/He_TACSham28d_RNA_XistBxC_-+_XX|XistBxC -+, CAST = Xi"
  "He_TAC|$ROOT/F1_TAC_Sarah/He_TAC28d_RNA_XistBxC_-+_XX|XistBxC -+, CAST = Xi"
)

# Batch: the array task picks the condition. Interactive: set IDX in the shell BEFORE
# running this block (IDX=3), or pass it per run (IDX=3 bash slurm/rnaseq_allelic_tracks.slurm).
IDX=${IDX:-${SLURM_ARRAY_TASK_ID:-1}}

# IDX may be a number (1..4) or a label from the table above, e.g. IDX=He_TAC
if [[ ! $IDX =~ ^[0-9]+$ ]]; then
  for i in "${!CONDITIONS[@]}"; do
    [[ ${CONDITIONS[$i]%%|*} == "$IDX" ]] && { IDX=$((i+1)); break; }
  done
fi
if [[ ! $IDX =~ ^[0-9]+$ ]] || (( IDX < 1 || IDX > ${#CONDITIONS[@]} )); then
  echo "IDX '$IDX' is not 1..${#CONDITIONS[@]} or one of: $(printf '%s ' "${CONDITIONS[@]%%|*}")" >&2
fi
IFS='|' read -r COND CDIR GENO <<< "${CONDITIONS[$((IDX-1))]}"

BAM=$CDIR/pooled.sorted.bam                             # strands merged: allele-only tracks
REF=$ROOT/bwa-index/GRCm39.primary_assembly.genome.fa
SNPBED=$ROOT/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed
VCF=$ROOT/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.vcf.gz   # the same VCF the ChIP script builds
OUT=$ROOT/heart_bulk_allelic_igv                        # shared by all four conditions
TMP=${SCRATCH:-$CDIR}/whatshap_${COND}

NPROC=${SLURM_CPUS_PER_TASK:-8}
JOBS=$(( NPROC < 8 ? NPROC : 8 ))
THRD_PER_JOB=$(( NPROC / JOBS )); THRD_PER_JOB=$(( THRD_PER_JOB > 0 ? THRD_PER_JOB : 1 ))

mkdir -p "$OUT" "$TMP"

conda activate RNAseq
for tool in whatshap samtools bgzip tabix bamCoverage; do
  command -v "$tool" >/dev/null || echo "MISSING TOOL: $tool is not in the RNAseq env" >&2
done

# He_TACSham28d has pooled_fwd/pooled_rev but no combined BAM: separate_BAM_strand.pl
# sends every read to exactly one of the two, so their union is the full library.
if [[ ! -f $BAM ]]; then
  echo "no pooled.sorted.bam in $CDIR - merging the strand BAMs"
  samtools merge -@ "$NPROC" -o "$TMP/pooled.bam" \
    "$CDIR/pooled_fwd.sorted.bam" "$CDIR/pooled_rev.sorted.bam"
  samtools sort -@ "$NPROC" -m 1G -o "$BAM" "$TMP/pooled.bam"
  rm -f "$TMP/pooled.bam"
fi
[[ -f ${BAM}.bai ]] || samtools index -@ "$NPROC" "$BAM"

echo "condition $IDX/${#CONDITIONS[@]}: $COND   [$GENO]"
echo "reads   $BAM"
echo "cores   $NPROC -> $JOBS haplotag workers x $THRD_PER_JOB output threads"
echo "scratch $TMP"
echo "output  $OUT"


# ==================================================================================
# BLOCK 1  -  VCF + reference checks. The VCF is shared with the ChIP script.
# ==================================================================================
# See chipseq_allelic_tracks.slurm for why column 2 of the SNP bed is the 1-based POS
# and why GT 0|1 makes haplotype 1 = B6 and haplotype 2 = CAST.
if [[ -f ${VCF}.tbi ]]; then
  echo "VCF already built: $VCF"
else
  echo "building $VCF"
  {
    echo '##fileformat=VCFv4.2'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    awk 'BEGIN{OFS="\t"} $4 ~ /^[ACGT][ACGT]$/ {
           print $1, $2, ".", substr($4,1,1), substr($4,2,1), "50.0", ".", ".", "GT", "0|1"
         }' "$SNPBED"
  } | bgzip -c > "$VCF"
  tabix -p vcf "$VCF"
fi

[[ -f ${REF}.fai ]] || samtools faidx "$REF"

# These BAMs came from STAR against gencode M36 / GRCm39, not from bwa, so confirm the
# contig lengths agree with the fasta whatshap will realign against before trusting a
# single tag. A silent mismatch here would mis-assign every read.
mismatch=$(join -j1 \
  <(samtools view -H "$BAM" |
    awk '$1=="@SQ"{ sn=""; ln="";
                    for(i=2;i<=NF;i++){ if($i~/^SN:/) sn=substr($i,4); else if($i~/^LN:/) ln=substr($i,4) }
                    if(sn!="") print sn, ln }' | sort -k1,1) \
  <(cut -f1,2 "${REF}.fai" | sort -k1,1) \
  | awk '$2 != $3 { print $1, "BAM="$2, "fasta="$3 }' | head -5)
if [[ -n $mismatch ]]; then
  echo "CONTIG LENGTH MISMATCH between BAM header and ${REF}.fai:" >&2
  echo "$mismatch" >&2
else
  echo "BAM header and reference agree on contig lengths"
fi


# ==================================================================================
# BLOCK 2  -  haplotag, one chromosome per worker.  ~45-90 min.
# ==================================================================================
# Paired-end: whatshap assigns a fragment by query name, so a SNP seen by either mate
# tags both. Spliced alignments are fine - reads are compared to the reference only
# where they actually align, so an intron simply contributes no SNPs.
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

{ for chrom in "${CHROMS[@]}"; do
    if [[ -s $TMP/${chrom}.haplotag.tsv ]]; then
      echo "  skip $chrom (list already exists)" >&2
    else
      echo "$chrom"
    fi
  done
} | xargs -n1 -P "$JOBS" -I{} bash -c 'tag_chrom "$@"' _ {}

{ head -n1 "$TMP/${CHROMS[0]}.haplotag.tsv"
  cat "$TMP"/*.haplotag.tsv | grep -v '^#'
} > "$TMP/haplotag_list.tsv"
echo "tagged alignments: H1 $(awk -F'\t' '$2=="H1"' "$TMP/haplotag_list.tsv" | wc -l), \
H2 $(awk -F'\t' '$2=="H2"' "$TMP/haplotag_list.tsv" | wc -l)"


# ==================================================================================
# BLOCK 3  -  split into the two alleles, sort, index, count.  ~40 min.
# ==================================================================================
whatshap split \
  --output-h1 "$TMP/${COND}_B6.bam" \
  --output-h2 "$TMP/${COND}_CAST.bam" \
  "$BAM" "$TMP/haplotag_list.tsv"

for allele in B6 CAST; do
  samtools sort -@ "$NPROC" -m 1G -o "$OUT/${COND}_${allele}.bam" "$TMP/${COND}_${allele}.bam"
  samtools index -@ "$NPROC" "$OUT/${COND}_${allele}.bam"
  rm -f "$TMP/${COND}_${allele}.bam"
done

N_B6=$(samtools view -c -@ "$NPROC" "$OUT/${COND}_B6.bam")
N_CAST=$(samtools view -c -@ "$NPROC" "$OUT/${COND}_CAST.bam")
SF=$(awk -v a="$N_B6" -v b="$N_CAST" 'BEGIN{ printf "%.10f", 1000000/(a+b) }')
echo "$COND assigned reads: B6 $N_B6, CAST $N_CAST"
awk -v a="$N_B6" -v b="$N_CAST" -v c="$COND" \
  'BEGIN{ printf "%s genome-wide B6 fraction: %.4f\n", c, a/(a+b) }'
echo "shared scale factor: $SF"
# Each condition gets its own denominator, so the two alleles are comparable within a
# condition and the conditions are comparable to each other in CPM. STAR ran with
# --outFilterMultimapNmax 1, so there are no multi-mappers to filter out here.


# ==================================================================================
# BLOCK 4  -  bigWigs.  ~20 min.  (needs $SF from BLOCK 3)
# ==================================================================================
# No --extendReads: these are spliced alignments, and extending a read across an intron
# would paint coverage over sequence that was never transcribed into that fragment.
bamCoverage -b "$OUT/${COND}_B6.bam"   -o "$OUT/${COND}_B6.bw" \
  --binSize "$BINSIZE" --scaleFactor "$SF" -p "$NPROC"
bamCoverage -b "$OUT/${COND}_CAST.bam" -o "$OUT/${COND}_CAST.bw" \
  --binSize "$BINSIZE" --scaleFactor "$SF" -p "$NPROC"


# ==================================================================================
# BLOCK 5  -  IGV session, one overlay row per condition. Run once, at the end.
# ==================================================================================
# Written to a temp file and moved into place, so parallel array tasks cannot interleave.
{
  cat <<XML
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" locus="chr5:148922383-148998282" version="8">
  <Resources>
XML
  for entry in "${CONDITIONS[@]}"; do
    IFS='|' read -r c _ _ <<< "$entry"
    printf '    <Resource path="%s_B6.bw"/>\n    <Resource path="%s_CAST.bw"/>\n' "$c" "$c"
  done
  echo '  </Resources>'
  echo '  <Panel name="DataPanel">'
  for entry in "${CONDITIONS[@]}"; do
    IFS='|' read -r c _ g <<< "$entry"
    cat <<XML
    <Track clazz="org.broad.igv.track.MergedTracks" id="${c}_allelic" attributeKey="${c}"
           name="${c}  [${g}]  B6 blue / CAST red" alpha="0.5" autoScale="false"
           renderer="BAR_CHART" height="70" visible="true">
      <DataRange minimum="0.0" baseline="0.0" maximum="${MAXV}.0" type="LINEAR"/>
      <Track clazz="org.broad.igv.track.DataSourceTrack" id="${c}_B6.bw" attributeKey="${c}_B6.bw"
             name="${c} B6" color="0,90,181" autoScale="false"
             renderer="BAR_CHART" windowFunction="mean" visible="true">
        <DataRange minimum="0.0" baseline="0.0" maximum="${MAXV}.0" type="LINEAR"/>
      </Track>
      <Track clazz="org.broad.igv.track.DataSourceTrack" id="${c}_CAST.bw" attributeKey="${c}_CAST.bw"
             name="${c} CAST" color="220,50,32" autoScale="false"
             renderer="BAR_CHART" windowFunction="mean" visible="true">
        <DataRange minimum="0.0" baseline="0.0" maximum="${MAXV}.0" type="LINEAR"/>
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
} > "$OUT/.heart_bulk_allelic_session.xml.tmp"
mv "$OUT/.heart_bulk_allelic_session.xml.tmp" "$OUT/heart_bulk_allelic_session.xml"

echo "Wrote:"
ls -lh "$OUT"
echo "scratch left in place for reruns (delete when happy): $TMP"
date
