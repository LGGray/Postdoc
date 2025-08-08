#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J RNAseq
#SBATCH -D  /dss/dssfs02/lwp-dss-0001/pr23fa/pr23fa-dss-0000/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --array=1-12

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs02/lwp-dss-0001/pr23fa/pr23fa-dss-0000/go93qiw2/miRNA/cardiac_miRNA

# ### Trim adaptor sequenecs from reads - TruSeq Small RNA 3'adaptor ###
# ### Ensure reads have a minimum read length as 18 bp and maximum as 30 bp
# for f in *_R1.fastq.gz; do
#     base="${f%.fastq.gz}"
#     cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 30 -o "${base}_trimmed.fastq.gz" "$f"
# done

# ### Align to mirBase Mouse miRNA sequences ###
# ### Keep only mouse sequences and replace U with T ###
# awk '/^>/ {keep = ($0 ~ /Mus musculus/)} keep' mature.fa > mouse_mature.fa
# sed 's/U/T/g' mouse_mature.fa > mouse_mature_DNA.fa
# bowtie-build mouse_mature_DNA.fa mouse_mature_DNA

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

bowtie \
  -n 0 -l 30 --norc --best --strata -m 1 -S \
  -x mouse_mature_DNA \
  "${sample}_R1_trimmed.fastq.gz" \
  --un "${sample}_unmapped_to_miRBase.fastq" \
  > "${sample}_miRBase.sam"

samtools view -bS "${sample}_miRBase.sam" | samtools sort -o "${sample}_miRBase.sorted.bam"
samtools index "${sample}_miRBase.sorted.bam"
samtools idxstats "${sample}_miRBase.sorted.bam" | cut -f1,3 > "${sample}_miRBase.counts.tsv"

# ### Take unaligned reads and align to the reference genome ###
# bowtie-build ../../genomeDir/GRCm38.primary_assembly.genome.fa GRCm38.primary_assembly.genome

bowtie \
  -n 1 -l 30 --norc --best --strata -m 1 -S \
  -x GRCm38.primary_assembly.genome \
  "${sample}_unmapped_to_miRBase.fastq" \
  --un "${sample}_unmapped_to_genome.fastq" \
  > "${sample}_genome.sam"

samtools view -bS "${sample}_genome.sam" | samtools sort -o "${sample}_genome.sorted.bam"
samtools index "${sample}_genome.sorted.bam"

# # Convert mmu.gff to bed file
# awk -v OFS="\t" '
#   $0 !~ /^#/ && $3=="miRNA" {
#     # parse the attributes column to get Name=...
#     name="";
#     n=split($9,a,";");
#     for(i=1;i<=n;i++){
#       split(a[i],kv,"=");
#       if(kv[1]=="Name"){ name=kv[2]; break }
#     }
#     if(name==""){        # fallback to ID if Name missing
#       for(i=1;i<=n;i++){ split(a[i],kv,"="); if(kv[1]=="ID"){ name=kv[2]; break } }
#     }
#     print $1, $4-1, $5, name, 0, $7
#   }
# ' ../mmu.gff3 \
# | sort -k1,1 -k2,2n > miRBase_mmu_mature.GRCm38.bed

# Tag reads with miRNA names when they overlap miRBase regions
bedtools tag -s -tag XT -i "${sample}_genome.sorted.bam" \
  -files miRBase_mmu_mature.GRCm38.bed \
  -names \
> "${sample}_genome.tagged.bam"
samtools index "${sample}_genome.tagged.bam"

samtools view "${sample}_genome.tagged.bam" | \
  awk '{
    tag="";
    for(i=12;i<=NF;i++){
      if($i ~ /^XT:Z:/){
        split($i,a,":");
        tag=a[3];
        counts[tag]++;
      }
    }
  }
  END{for(k in counts) print k"\t"counts[k]}
' > "${sample}_genome_miRNA.counts.tsv"

for f in *_miRBase.counts.tsv; do
  sample="${f%_miRBase.counts.tsv}"
  miR="${sample}_miRBase.counts.tsv"
  gen="${sample}_genome_miRNA.counts.tsv"

  # Make sure missing genome file becomes empty 0s
  [ -f "$gen" ] || touch "$gen"

  # Collapse duplicates, drop the '*' row, and sort
  awk -F'\t' '$1!="*"{c[$1]+=$2} END{for(k in c) print k"\t"c[k]}' "$miR" \
    | sort -k1,1 > "${sample}.miR.tmp"
  awk -F'\t' '$1!="*"{c[$1]+=$2} END{for(k in c) print k"\t"c[k]}' "$gen" \
    | sort -k1,1 > "${sample}.gen.tmp"

  # Outer-join and sum counts
  join -t $'\t' -a1 -a2 -e 0 -o 0,1.2,2.2 \
    "${sample}.miR.tmp" "${sample}.gen.tmp" \
  | awk -F'\t' -v OFS="\t" -v S="$sample" 'BEGIN{print "miRNA", S} {print $1, $2+$3}' \
  > "${sample}.final.tsv"

  rm -f "${sample}.miR.tmp" "${sample}.gen.tmp"
done
