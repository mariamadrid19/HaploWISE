#!/bin/bash -l 
#SBATCH --cluster=wice
#SBATCH --job-name call_snps
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o call_snps.%j.out
#SBATCH --array=1-11

set -euo pipefail

# Index starts at 0
ID=$((SLURM_ARRAY_TASK_ID - 1))

# Load required modules
module load BCFtools/1.12-GCC-10.3.0 
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate variant_tools

# Reference files
REFNAME="Pchal_Bar_LW"
REF="/scratch/leuven/357/vsc35707/popgen/genomes/Pchalceus_LW_final.fasta"
FAI="${REF}.fai"

# Output directories
VCF_DIR="vcfs"
OUT_DIR="final-calls"
mkdir -p "$VCF_DIR"
mkdir -p "$OUT_DIR"

echo "================="

# Sample IDs (all samples, 123)
samples=(\
GC129388 GC129389 GC129390 GC129391 GC129392 \
GC129393 GC129394 GC129395 GC129396 GC129397 \
GC129398 GC129399 GC129400 GC129401 GC129402 \
GC129403 GC129404 GC129405 GC129406 GC129407 \
GC129408 GC129409 GC129410 GC129411 GC129412 \
GC129413 GC129414 GC129415 GC129416 GC129417 \
GC129418 GC129419 GC129420 GC129421 GC129422 \
GC129423 GC129424 GC129425 GC129426 GC129427 \
GC129428 GC129429 GC129430 GC129431 GC129432 \
GC129433 GC129434 GC129435 GC129437 GC129438 \
GC129439 GC129440 GC136078 GC136079 GC136080 \
GC136081 GC136082 GC136083 GC136084 GC136085 \
GC136086 GC136087 GC136088 GC136089 GC136090 \
GC136091 GC136092 GC136093 GC136094 GC136095 \
GC136096 GC136097 GC136098 GC136099 GC136100 \
GC136101 GC136102 GC136103 GC136104 GC136105 \
GC136106 GC136107 GC136108 GC136109 GC136110 \
GC136111 GC136112 GC136113 GC136114 GC136115 \
GC136116 GC136117 GC136118 GC136119 GC136120 \
GC136121 GC136122 GC136123 GC136124 GC136125 \
GC136126 GC136127 GC136128 SW_01 SW_02 SW_03 SW_04 SW_05 \
SW_06 SW_07 SW_08 SW_09 SW_10 \
LW_01 LW_02 LW_03 LW_04 LW_05 \
LW_06 LW_07 LW_08 LW_09 LW_10)


CHR=(CHR1 CHR2 CHR3 CHR4 CHR5 CHR6 CHR7 CHR8 CHR9 CHR10 CHR11)
names=(1 2 3 4 5 6 7 8 9 10 11)

BAM_DIR="/scratch/leuven/357/vsc35707/popgen/bams"

ALL_LIST=()

for SAMPLE in "${samples[@]}"; do
    BAM="${BAM_DIR}/${SAMPLE}.${REFNAME}.filtered.sorted.dedup.bam"

    if [[ ! -f "$BAM" ]]; then
        echo "ERROR: Missing BAM: $BAM" >&2
        exit 1
    fi

    ALL_LIST+=("$BAM")
done

# Step 1: Call variants with FORMAT annotations
bcftools mpileup \
    -Ou \
    --threads 36 \
    --fasta-ref "$REF" \
    --regions "${CHR[ID]}" \
    --annotate FORMAT/DP,FORMAT/AD \
    "${ALL_LIST[@]}" | \
bcftools call \
    --multiallelic-caller \
    -Oz \
    -f GQ \
    -o "$VCF_DIR/$REFNAME.chr_${names[ID]}.vcf.gz"

# Step 2: Filter with vcftools and index
vcftools --gzvcf $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").vcf.gz --recode --remove-indels --stdout | bgzip > $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").filtered.vcf.gz

bcftools index -t "$VCF_DIR/$REFNAME.chr_${names[ID]}.filtered.vcf.gz"

# Step 3: Phase with beagle without imputation
module load Java/21.0.8

java -Xmx32g -jar /data/leuven/357/vsc35707/beagle.27Feb25.75f.jar \
  gt=$VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").filtered.vcf.gz \
  out=$VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").phased \
  nthreads=36

# Step 4: Parse VCF with custom script
/data/leuven/357/vsc35707/miniconda3/bin/python3.11 parseVCF.py \
  --gtf flag=GQ min=30 gtTypes=Het \
  --gtf flag=GQ min=30 gtTypes=HomAlt \
  --gtf flag=DP min=10 \
  --skipIndels \
  -i $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").phased.vcf.gz \
  | gzip > $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").calls.gz

echo "Calls file created: $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").calls.gz" 

zcat "$VCF_DIR/$REFNAME.chr_${names[ID]}.calls.gz" | \
awk -v OFS="\t" 'NR==1 {
    for (i = 1; i <= NF; i++) {
        gsub(/^.*\//, "", $i);
        gsub(/\.Pchal_Bar_LW\.filtered\.sorted\.dedup\.bam$/, "", $i);
    }
}
{ print }
' | gzip -c > "$OUT_DIR/$REFNAME.chr_${names[ID]}.H.calls.gz"

echo "Header fixed and saved: $OUT_DIR/$REFNAME.chr_${names[ID]}.H.calls.gz"
