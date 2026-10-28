#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --time=60:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o /scratch/leuven/357/vsc35707/popgen/scripts/logs/map_all.%A_%a.out
#SBATCH --array=1-123%20

set -euo pipefail

# =========================================================
# 1. Settings
# =========================================================

ID=$((SLURM_ARRAY_TASK_ID - 1))
THREADS="${SLURM_CPUS_PER_TASK:-36}"

# Activate Conda environment
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate variant_tools

# Load Picard
module load picard/3.4.0-Java-17

# Directories
BASE_DIR="/scratch/leuven/357/vsc35707/popgen"
REF="${BASE_DIR}/genomes/Pchalceus_LW_final.fasta"
REFNAME="Pchal_Bar_LW"

READS_DIR="${BASE_DIR}/reads"
BAM_DIR="${BASE_DIR}/bams"

mkdir -p "$BAM_DIR"

# =========================================================
# 2. Sample IDs
# =========================================================

samples=(
GC129388 GC129389 GC129390 GC129391 GC129392
GC129393 GC129394 GC129395 GC129396 GC129397
GC129398 GC129399 GC129400 GC129401 GC129402
GC129403 GC129404 GC129405 GC129406 GC129407
GC129408 GC129409 GC129410 GC129411 GC129412
GC129413 GC129414 GC129415 GC129416 GC129417
GC129418 GC129419 GC129420 GC129421 GC129422
GC129423 GC129424 GC129425 GC129426 GC129427
GC129428 GC129429 GC129430 GC129431 GC129432
GC129433 GC129434 GC129435 GC129437 GC129438
GC129439 GC129440 GC136078 GC136079 GC136080
GC136081 GC136082 GC136083 GC136084 GC136085
GC136086 GC136087 GC136088 GC136089 GC136090
GC136091 GC136092 GC136093 GC136094 GC136095
GC136096 GC136097 GC136098 GC136099 GC136100
GC136101 GC136102 GC136103 GC136104 GC136105
GC136106 GC136107 GC136108 GC136109 GC136110
GC136111 GC136112 GC136113 GC136114 GC136115
GC136116 GC136117 GC136118 GC136119 GC136120
GC136121 GC136122 GC136123 GC136124 GC136125
GC136126 GC136127 GC136128
SW_01 SW_02 SW_03 SW_04 SW_05
SW_06 SW_07 SW_08 SW_09 SW_10
LW_01 LW_02 LW_03 LW_04 LW_05
LW_06 LW_07 LW_08 LW_09 LW_10
)

if [[ ${#samples[@]} -ne 123 ]]; then
    echo "ERROR: Expected 123 samples, found ${#samples[@]}"
    exit 1
fi

SAMPLE="${samples[$ID]}"

# =========================================================
# 3. Input and output files
# =========================================================

FILE1="${READS_DIR}/${SAMPLE}_R1.fastq.gz"
FILE2="${READS_DIR}/${SAMPLE}_R2.fastq.gz"

RAW_BAM="${BAM_DIR}/${SAMPLE}.${REFNAME}.bam"
FILTERED_BAM="${BAM_DIR}/${SAMPLE}.${REFNAME}.filtered.bam"
SORTED_BAM="${BAM_DIR}/${SAMPLE}.${REFNAME}.filtered.sorted.bam"
DEDUP_BAM="${BAM_DIR}/${SAMPLE}.${REFNAME}.filtered.sorted.dedup.bam"

METRICS="${BAM_DIR}/${SAMPLE}.${REFNAME}.dup_metrics.txt"

TMPDIR_SAMPLE="${BAM_DIR}/temp_${SAMPLE}"

mkdir -p "$TMPDIR_SAMPLE"

echo "======================================"
echo "Sample: $SAMPLE"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Threads: $THREADS"
echo "======================================"

# Check input files
if [[ ! -s "$FILE1" || ! -s "$FILE2" ]]; then
    echo "ERROR: Missing FASTQ files for $SAMPLE" >&2
    exit 1
fi

if [[ ! -s "$REF" ]]; then
    echo "ERROR: Reference genome not found: $REF" >&2
    exit 1
fi

# Avoid overwriting completed results
if [[ -s "$DEDUP_BAM" && -s "${DEDUP_BAM}.bai" ]]; then
    if samtools quickcheck "$DEDUP_BAM"; then
        echo "Final BAM already exists and passed quickcheck."
        echo "Skipping $SAMPLE"
        exit 0
    fi
fi

# =========================================================
# 4. Map reads with BWA MEM
# =========================================================

echo "Mapping reads: $SAMPLE"

# Add read groups directly during alignment.
# Each sample receives a unique read-group ID and sample name.

RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:${SAMPLE}"

bwa mem \
    -t "$THREADS" \
    -M \
    -R "$RG" \
    "$REF" \
    "$FILE1" \
    "$FILE2" | \
samtools view \
    -@ "$THREADS" \
    -b \
    -o "$RAW_BAM" -

samtools quickcheck -v "$RAW_BAM"

# =========================================================
# 5. Filter reads
# =========================================================

echo "Filtering reads: $SAMPLE"

# Keep properly paired reads with MAPQ >= 20.

samtools view \
    -@ "$THREADS" \
    -f 0x02 \
    -q 20 \
    -b \
    -o "$FILTERED_BAM" \
    "$RAW_BAM"

samtools quickcheck -v "$FILTERED_BAM"

# Remove raw BAM after successful filtering
if [[ $(stat -c %s "$FILTERED_BAM") -gt 512000 ]]; then
    rm "$RAW_BAM"
fi

# =========================================================
# 6. Sort BAM
# =========================================================

echo "Sorting reads: $SAMPLE"

samtools sort \
    -@ "$THREADS" \
    -T "${TMPDIR_SAMPLE}/sort" \
    -o "$SORTED_BAM" \
    "$FILTERED_BAM"

samtools quickcheck -v "$SORTED_BAM"

# Remove filtered BAM after successful sorting
if [[ $(stat -c %s "$SORTED_BAM") -gt 512000 ]]; then
    rm "$FILTERED_BAM"
fi

# =========================================================
# 7. Remove PCR duplicates
# =========================================================

echo "Removing PCR duplicates: $SAMPLE"

java -Xmx8g \
    -Djava.io.tmpdir="$TMPDIR_SAMPLE" \
    -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
    I="$SORTED_BAM" \
    O="$DEDUP_BAM" \
    REMOVE_DUPLICATES=true \
    M="$METRICS" \
    ASSUME_SORTED=true \
    TMP_DIR="$TMPDIR_SAMPLE"

samtools quickcheck -v "$DEDUP_BAM"

# =========================================================
# 8. Index final BAM
# =========================================================

echo "Indexing BAM: $SAMPLE"

samtools index \
    -@ "$THREADS" \
    "$DEDUP_BAM"

# =========================================================
# 9. Validate final output and clean up
# =========================================================

if [[ -s "$DEDUP_BAM" && -s "${DEDUP_BAM}.bai" ]]; then

    samtools quickcheck -v "$DEDUP_BAM"

    echo "Final BAM and index successfully created."

    # Remove intermediate BAMs
    rm -f "$RAW_BAM"
    rm -f "$FILTERED_BAM"
    rm -f "$SORTED_BAM"

    # Remove temporary directory
    rm -rf -- "$TMPDIR_SAMPLE"

else
    echo "ERROR: Final BAM or index missing for $SAMPLE" >&2
    exit 1
fi

echo "======================================"
echo "Finished sample: $SAMPLE"
echo "Final BAM: $DEDUP_BAM"
echo "======================================"
