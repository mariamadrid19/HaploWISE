#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=call_snps
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=72G
#SBATCH --time=72:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o /scratch/leuven/357/vsc35707/popgen/scripts/logs/call_snps.%A_%a.out
#SBATCH --array=1-11%3

set -euo pipefail

# =========================================================
# 1. SETTINGS
# =========================================================

BASE_DIR="/scratch/leuven/357/vsc35707/popgen"

cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-12}"
ID=$((SLURM_ARRAY_TASK_ID - 1))

# Environment
module load BCFtools/1.12-GCC-10.3.0

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate variant_tools

module load Java/21.0.8

# Reference
REFNAME="Pchal_Bar_LW"
REF="${BASE_DIR}/genomes/Pchalceus_LW_final.fasta"

# Directories
BAM_DIR="${BASE_DIR}/bams"
VCF_DIR="${BASE_DIR}/vcfs"
OUT_DIR="${BASE_DIR}/final-calls"

mkdir -p "$VCF_DIR" "$OUT_DIR"

# Chromosomes
CHR=(CHR1 CHR2 CHR3 CHR4 CHR5 CHR6 CHR7 CHR8 CHR9 CHR10 CHR11)

CHROM="${CHR[$ID]}"
CHR_NUM=$((ID + 1))

PREFIX="${VCF_DIR}/${REFNAME}.chr_${CHR_NUM}"

FINAL="${OUT_DIR}/${REFNAME}.chr_${CHR_NUM}.H.calls.gz"

if [[ -s "$FINAL" ]] && gzip -t "$FINAL" 2>/dev/null; then
    echo "$CHROM already has a valid compressed H calls file."
    echo "Skipping to avoid overwriting existing results."
    exit 0
fi

# Filtering thresholds
MIN_DP=10
MIN_GQ=30

echo "======================================"
echo "Chromosome: $CHROM"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Threads: $THREADS"
echo "Minimum DP: $MIN_DP"
echo "Minimum GQ (Het/HomAlt): $MIN_GQ"
echo "======================================"

# =========================================================
# 2. SAMPLE IDS
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

# =========================================================
# 3. CHECK INPUTS AND SOFTWARE
# =========================================================

echo "Checking input BAM files..."

# Python dependency for filtering and restoration
if ! python3.11 -c "import pysam" 2>/dev/null; then
    echo "ERROR: pysam is not available in Python 3.11." >&2
    echo "Install pysam in the variant_tools environment first." >&2
    exit 1
fi

if [[ ! -s "$REF" || ! -s "${REF}.fai" ]]; then
    echo "ERROR: Reference FASTA or FASTA index missing." >&2
    exit 1
fi

if ! awk -v chr="$CHROM" '
    $1 == chr {found=1}
    END {exit !found}
' "${REF}.fai"; then
    echo "ERROR: $CHROM missing from reference index." >&2
    exit 1
fi

ALL_LIST=()

for SAMPLE in "${samples[@]}"; do

    BAM="${BAM_DIR}/${SAMPLE}.${REFNAME}.filtered.sorted.dedup.bam"

    if [[ ! -s "$BAM" ]]; then
        echo "ERROR: Missing BAM: $BAM" >&2
        exit 1
    fi

    if [[ ! -s "${BAM}.bai" ]]; then
        echo "ERROR: Missing BAM index: ${BAM}.bai" >&2
        exit 1
    fi

    if ! samtools quickcheck "$BAM"; then
        echo "ERROR: Invalid BAM: $BAM" >&2
        exit 1
    fi

    ALL_LIST+=("$BAM")

done

echo "BAM files found: ${#ALL_LIST[@]}"

if [[ "${#ALL_LIST[@]}" -ne 123 ]]; then
    echo "ERROR: Expected 123 BAM files." >&2
    exit 1
fi

# =========================================================
# 4. CALL VARIANTS
# =========================================================

echo "Calling variants on $CHROM..."

bcftools mpileup \
    -Ou \
    --threads "$THREADS" \
    --fasta-ref "$REF" \
    --regions "$CHROM" \
    --annotate FORMAT/DP,FORMAT/AD \
    "${ALL_LIST[@]}" | \
bcftools call \
    --multiallelic-caller \
    -Oz \
    -f GQ \
    -o "${PREFIX}.vcf.gz"

bcftools index -t "${PREFIX}.vcf.gz"

# Verify sample count
N_SAMPLES=$(bcftools query -l "${PREFIX}.vcf.gz" | wc -l)

if [[ "$N_SAMPLES" -ne "${#samples[@]}" ]]; then
    echo "ERROR: VCF contains $N_SAMPLES samples; expected ${#samples[@]}." >&2
    exit 1
fi

echo "Variant calling completed."

# =========================================================
# 5. REMOVE INDELS
# =========================================================

echo "Removing indels..."

bcftools view \
    -v snps \
    -Oz \
    -o "${PREFIX}.snps.vcf.gz" \
    "${PREFIX}.vcf.gz"

bcftools index -t "${PREFIX}.snps.vcf.gz"

echo "SNP-only VCF created."

# =========================================================
# 6. CREATE PYTHON HELPER
# =========================================================

# This helper performs two operations:
#
# filter:
#   Mask low-quality genotypes BEFORE phasing.
#
# restore:
#   Restore original missingness AFTER phasing,
#   while preserving Beagle's phase assignments.
#
# Each chromosome job creates its own helper in a
# task-specific temporary directory.

TMPDIR_SAMPLE="${BASE_DIR}/tmp_call_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$TMPDIR_SAMPLE"
trap 'status=$?; if [[ $status -eq 0 ]]; then rm -rf -- "$TMPDIR_SAMPLE"; else echo "Temporary files retained: $TMPDIR_SAMPLE" >&2; fi' EXIT

HELPER="${TMPDIR_SAMPLE}/process_genotypes.py"

cat > "$HELPER" <<'PYTHON'
#!/usr/bin/env python3
import argparse
from itertools import zip_longest
import pysam


def key(rec):
    return (rec.contig, rec.pos, rec.ref, tuple(rec.alts or ()))


def missing(gt):
    return gt is None or any(a is None for a in gt)


def filter_calls(src, dest, min_dp, min_gq):
    with pysam.VariantFile(src) as inp:
        out = pysam.VariantFile(dest, 'wz', header=inp.header)
        try:
            samples = list(inp.header.samples)
            counts = dict(total=0, original_missing=0, low_dp=0, low_gq=0)
            for rec in inp:
                for sample in samples:
                    call = rec.samples[sample]
                    gt = call.get('GT')
                    counts['total'] += 1
                    if missing(gt):
                        counts['original_missing'] += 1
                        continue
                    dp = call.get('DP')
                    if dp is None or dp < min_dp:
                        call['GT'] = (None,) * len(gt)
                        call.phased = False
                        counts['low_dp'] += 1
                    elif any(a != 0 for a in gt):
                        gq = call.get('GQ')
                        if gq is None or gq < min_gq:
                            call['GT'] = (None,) * len(gt)
                            call.phased = False
                            counts['low_gq'] += 1
                out.write(rec)
        finally:
            out.close()
    print('Filtering:', counts, flush=True)


def restore(src, beagle, dest):

    # The Beagle header has already been repaired
    # with bcftools reheader before this function.

    with pysam.VariantFile(src) as original, \
         pysam.VariantFile(beagle) as phased:

        samples = list(original.header.samples)

        if samples != list(phased.header.samples):
            raise ValueError(
                "Filtered and Beagle sample names/order differ"
            )

        if len(samples) != 123:
            raise ValueError(
                f"Expected 123 samples, found {len(samples)}"
            )

        header = phased.header

        restored = 0
        already_missing = 0
        retained = 0
        variants = 0

        with pysam.VariantFile(
            dest, 'wz', header=header
        ) as out:

            for orig, phase in zip_longest(original, phased):

                if orig is None or phase is None:
                    raise ValueError(
                        "Filtered and Beagle record counts differ"
                    )

                if key(orig) != key(phase):
                    raise ValueError(
                        f"Variant mismatch: {key(orig)} vs {key(phase)}"
                    )

                variants += 1

                for sample in samples:

                    gt = orig.samples[sample].get('GT')
                    call = phase.samples[sample]
                    phased_gt = call.get('GT')

                    if missing(gt):

                        if missing(phased_gt):
                            already_missing += 1
                        else:
                            restored += 1

                        ploidy = len(gt) if gt is not None else 2

                        call['GT'] = (None,) * ploidy
                        call.phased = False

                    else:

                        if missing(phased_gt):
                            raise ValueError(
                                f"Observed genotype missing at "
                                f"{key(orig)}, {sample}"
                            )

                        if sorted(gt) != sorted(phased_gt):
                            raise ValueError(
                                f"Observed genotype changed at "
                                f"{key(orig)}, {sample}"
                            )

                        retained += 1

                # Input records and output use the
                # same repaired Beagle header.
                out.write(phase)

                if variants % 100000 == 0:
                    print(
                        f"Restored {variants:,} variant records",
                        flush=True
                    )

        print(
            f"Restoration complete: {variants:,} variants, "
            f"{retained:,} observed retained, "
            f"{restored:,} Beagle-filled masked, "
            f"{already_missing:,} already missing",
            flush=True
        )

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['filter', 'restore'])
    p.add_argument('--input')
    p.add_argument('--filtered')
    p.add_argument('--phased')
    p.add_argument('--output', required=True)
    p.add_argument('--min-dp', type=int, default=10)
    p.add_argument('--min-gq', type=int, default=30)
    a = p.parse_args()
    if a.mode == 'filter':
        filter_calls(a.input, a.output, a.min_dp, a.min_gq)
    else:
        restore(a.filtered, a.phased, a.output)
PYTHON

python3.11 -m py_compile "$HELPER"

# =========================================================
# 7. FILTER GENOTYPES BEFORE PHASING
# =========================================================

echo "Filtering genotypes before phasing..."

python3.11 "$HELPER" filter \
    --input "${PREFIX}.snps.vcf.gz" \
    --output "${PREFIX}.filtered.vcf.gz" \
    --min-dp "$MIN_DP" \
    --min-gq "$MIN_GQ"

bcftools index -t "${PREFIX}.filtered.vcf.gz"

echo "Filtered VCF created:"
echo "${PREFIX}.filtered.vcf.gz"

# =========================================================
# 8. PHASE WITH BEAGLE
# =========================================================

echo "Phasing $CHROM with Beagle..."

java -Xmx64g \
    -jar /data/leuven/357/vsc35707/beagle.27Feb25.75f.jar \
    gt="${PREFIX}.filtered.vcf.gz" \
    out="${PREFIX}.beagle" \
    nthreads=12 \
    impute=false \
    window=5.0 \
    overlap=1.0

echo "Beagle phasing completed."

# =========================================================
# 9. REPAIR BEAGLE HEADER AND RESTORE MISSINGNESS
# =========================================================

echo "Repairing Beagle VCF header..."

# Beagle can omit ##contig definitions.
# Repair the header before pysam reads variant records.

FIXED_BEAGLE="${TMPDIR_SAMPLE}/beagle.fixed.vcf.gz"

bcftools reheader \
    -f "${REF}.fai" \
    -o "$FIXED_BEAGLE" \
    "${PREFIX}.beagle.vcf.gz"

# Check that the chromosome definition exists.

if ! bcftools view -h "$FIXED_BEAGLE" |
    grep -F "##contig=<ID=${CHROM}," > /dev/null; then

    echo "ERROR: Missing contig definition for $CHROM" >&2
    exit 1
fi

# Index repaired Beagle output.

bcftools index -f -t "$FIXED_BEAGLE"

# Validate variant counts before restoration.

N_FILTERED=$(bcftools index -n "${PREFIX}.filtered.vcf.gz")
N_BEAGLE=$(bcftools index -n "$FIXED_BEAGLE")

if [[ "$N_FILTERED" -ne "$N_BEAGLE" ]]; then
    echo "ERROR: Filtered and Beagle variant counts differ." >&2
    exit 1
fi

# Validate sample names and order.

bcftools query -l "${PREFIX}.filtered.vcf.gz" \
    > "${TMPDIR_SAMPLE}/filtered.samples"

bcftools query -l "$FIXED_BEAGLE" \
    > "${TMPDIR_SAMPLE}/beagle.samples"

if ! cmp -s \
    "${TMPDIR_SAMPLE}/filtered.samples" \
    "${TMPDIR_SAMPLE}/beagle.samples"; then

    echo "ERROR: Sample names/order differ." >&2
    exit 1
fi

echo "Beagle VCF validated."

# Restore original missingness.

echo "Restoring missing genotypes..."

TEMP_PHASED="${TMPDIR_SAMPLE}/phased.vcf.gz"

python3.11 "$HELPER" restore \
    --filtered "${PREFIX}.filtered.vcf.gz" \
    --phased "$FIXED_BEAGLE" \
    --output "$TEMP_PHASED"

# Index and validate final phased VCF.

bcftools index -f -t "$TEMP_PHASED"

N_FINAL=$(bcftools index -n "$TEMP_PHASED")
N_SAMPLES=$(bcftools query -l "$TEMP_PHASED" | wc -l)

if [[ "$N_FINAL" -ne "$N_FILTERED" ]]; then
    echo "ERROR: Final variant count mismatch." >&2
    exit 1
fi

if [[ "$N_SAMPLES" -ne 123 ]]; then
    echo "ERROR: Expected 123 samples." >&2
    exit 1
fi

mv -f "$TEMP_PHASED" "${PREFIX}.phased.vcf.gz"
mv -f "${TEMP_PHASED}.tbi" "${PREFIX}.phased.vcf.gz.tbi"

echo "Phased VCF created:"
echo "${PREFIX}.phased.vcf.gz"


# =========================================================
# 10. PARSE PHASED VCF
# =========================================================

echo "Parsing phased VCF..."

# DP and GQ filtering were already performed.
# Do not repeat these filters after Beagle.

python3.11 "${BASE_DIR}/parseVCF.py" \
    --skipIndels \
    -i "${PREFIX}.phased.vcf.gz" \
    | gzip -c > "${TMPDIR_SAMPLE}/calls.gz"
gzip -t "${TMPDIR_SAMPLE}/calls.gz"
mv -f "${TMPDIR_SAMPLE}/calls.gz" "${PREFIX}.calls.gz"

echo "Calls file created:"
echo "${PREFIX}.calls.gz"

# =========================================================
# 11. FIX SAMPLE NAMES IN CALLS HEADER
# =========================================================

echo "Fixing sample names..."

zcat "${PREFIX}.calls.gz" | \
awk -v OFS="\t" 'NR==1 {
    for (i = 1; i <= NF; i++) {
        gsub(/^.*\//, "", $i);
        gsub(/\.Pchal_Bar_LW\.filtered\.sorted\.dedup\.bam$/, "", $i);
    }
}
{ print }
' | gzip -c > "${TMPDIR_SAMPLE}/H.calls.gz"
gzip -t "${TMPDIR_SAMPLE}/H.calls.gz"

H_LINES=$(gzip -cd "${TMPDIR_SAMPLE}/H.calls.gz" |
    awk 'END {print NR}')

if [[ "$H_LINES" -le 1 ]]; then
    echo "ERROR: H calls file contains no variant data." >&2
    exit 1
fi

echo "H calls lines: $H_LINES"

mv -f "${TMPDIR_SAMPLE}/H.calls.gz" "${OUT_DIR}/${REFNAME}.chr_${CHR_NUM}.H.calls.gz"

echo "Final calls file:"
echo "${OUT_DIR}/${REFNAME}.chr_${CHR_NUM}.H.calls.gz"

# =========================================================
# 12. CLEAN UP TEMPORARY HELPER
# =========================================================

echo "======================================"
echo "FINISHED CHROMOSOME $CHROM"
echo "======================================"
