#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=call_snps
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=48G
#SBATCH --time=72:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o /scratch/leuven/357/vsc35707/popgen/scripts/logs/call_snps.%A_%a.out
#SBATCH --array=1-11

set -euo pipefail

# =========================================================
# 1. SETTINGS
# =========================================================

BASE_DIR="/scratch/leuven/357/vsc35707/popgen"

cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-36}"
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

HELPER="${TMPDIR_SAMPLE}/process_genotypes.py"

cat > "$HELPER" <<'PYTHON'

#!/usr/bin/env python3

import argparse
import sys
from itertools import zip_longest

import pysam


def variant_key(record):
    return (
        record.contig,
        record.pos,
        record.ref,
        tuple(record.alts or ())
    )


def is_missing(gt):
    return gt is None or any(a is None for a in gt)


def filter_genotypes(input_file, output_file, min_dp, min_gq):

    invcf = pysam.VariantFile(input_file)
    outvcf = pysam.VariantFile(
        output_file,
        "wz",
        header=invcf.header
    )

    samples = list(invcf.header.samples)

    total_genotypes = 0
    masked_dp = 0
    masked_gq = 0
    already_missing = 0

    for record in invcf:

        for sample_name in samples:

            call = record.samples[sample_name]
            gt = call.get("GT")

            total_genotypes += 1

            if is_missing(gt):
                already_missing += 1
                continue

            dp = call.get("DP")

            # DP filtering applies to every genotype.
            if dp is None or dp < min_dp:

                call["GT"] = tuple(None for _ in gt)
                call.phased = False

                masked_dp += 1
                continue

            # Apply GQ threshold to heterozygous
            # and homozygous-alternative genotypes.
            #
            # Homozygous-reference genotypes are not
            # filtered on GQ, matching the original
            # parseVCF.py filtering logic.

            is_hom_ref = all(a == 0 for a in gt)

            if not is_hom_ref:

                gq = call.get("GQ")

                if gq is None or gq < min_gq:

                    call["GT"] = tuple(None for _ in gt)
                    call.phased = False

                    masked_gq += 1

        outvcf.write(record)

    invcf.close()
    outvcf.close()

    print("Genotype filtering completed.")
    print(f"Total genotypes: {total_genotypes}")
    print(f"Already missing: {already_missing}")
    print(f"Masked due to DP: {masked_dp}")
    print(f"Masked due to GQ: {masked_gq}")
    print(
        "Total missing after filtering: "
        f"{already_missing + masked_dp + masked_gq}"
    )


def restore_missingness(filtered_file, phased_file, output_file):

    original = pysam.VariantFile(filtered_file)
    phased = pysam.VariantFile(phased_file)

    original_samples = list(original.header.samples)
    phased_samples = list(phased.header.samples)

    # Verify sample identities and order.
    if original_samples != phased_samples:
        raise RuntimeError(
            "Sample names or sample order differ between "
            "filtered VCF and Beagle VCF."
        )

    outvcf = pysam.VariantFile(
        output_file,
        "wz",
        header=phased.header
    )

    restored = 0
    retained = 0
    n_variants = 0

    for original_record, phased_record in zip_longest(
        original,
        phased
    ):

        if original_record is None or phased_record is None:
            raise RuntimeError(
                "Filtered and phased VCFs have different "
                "numbers of variant records."
            )

        if variant_key(original_record) != variant_key(phased_record):
            raise RuntimeError(
                "Variant mismatch between filtered and phased VCFs: "
                f"{variant_key(original_record)} vs "
                f"{variant_key(phased_record)}"
            )

        n_variants += 1

        for sample_name in original_samples:

            original_gt = original_record.samples[sample_name].get("GT")
            phased_call = phased_record.samples[sample_name]

            if is_missing(original_gt):

                phased_gt = phased_call.get("GT")

                # Restore missing genotype, including ploidy.
                phased_call["GT"] = tuple(
                    None for _ in original_gt
                )

                phased_call.phased = False

                if not is_missing(phased_gt):
                    restored += 1

            else:

                # Preserve Beagle's phased genotype.
                #
                # Also verify that Beagle did not change
                # the observed genotype's allele counts.

                phased_gt = phased_call.get("GT")

                if is_missing(phased_gt):
                    raise RuntimeError(
                        "Beagle produced a missing genotype "
                        "where the filtered input was observed: "
                        f"{sample_name}, "
                        f"{variant_key(original_record)}"
                    )

                if sorted(original_gt) != sorted(phased_gt):
                    raise RuntimeError(
                        "Beagle changed an observed genotype: "
                        f"{sample_name}, "
                        f"{variant_key(original_record)}"
                    )

                retained += 1

        outvcf.write(phased_record)

    original.close()
    phased.close()
    outvcf.close()

    print("Missingness restoration completed.")
    print(f"Variants checked: {n_variants}")
    print(f"Observed genotypes retained: {retained}")
    print(f"Beagle-filled genotypes masked again: {restored}")


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "mode",
        choices=["filter", "restore"]
    )

    parser.add_argument("--input")
    parser.add_argument("--filtered")
    parser.add_argument("--phased")

    parser.add_argument("--output", required=True)

    parser.add_argument("--min-dp", type=int, default=10)
    parser.add_argument("--min-gq", type=int, default=30)

    args = parser.parse_args()

    if args.mode == "filter":

        filter_genotypes(
            args.input,
            args.output,
            args.min_dp,
            args.min_gq
        )

    elif args.mode == "restore":

        restore_missingness(
            args.filtered,
            args.phased,
            args.output
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

PYTHON

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

java -Xmx32g \
    -jar /data/leuven/357/vsc35707/beagle.27Feb25.75f.jar \
    gt="${PREFIX}.filtered.vcf.gz" \
    out="${PREFIX}.beagle" \
    nthreads="$THREADS" \
    impute=false

echo "Beagle phasing completed."

# =========================================================
# 9. RESTORE ORIGINAL MISSINGNESS
# =========================================================

echo "Restoring missing genotypes..."

python3.11 "$HELPER" restore \
    --filtered "${PREFIX}.filtered.vcf.gz" \
    --phased "${PREFIX}.beagle.vcf.gz" \
    --output "${PREFIX}.phased.vcf.gz"

bcftools index -t "${PREFIX}.phased.vcf.gz"

echo "Final phased VCF created:"
echo "${PREFIX}.phased.vcf.gz"

# =========================================================
# 10. PARSE PHASED VCF
# =========================================================

echo "Parsing phased VCF..."

# DP and GQ filtering were already performed.
# Do not repeat these filters after Beagle.

/data/leuven/357/vsc35707/miniconda3/bin/python3.11 \
    parseVCF.py \
    --skipIndels \
    -i "${PREFIX}.phased.vcf.gz" \
    | gzip -c > "${PREFIX}.calls.gz"

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
' | gzip -c > "${OUT_DIR}/${REFNAME}.chr_${CHR_NUM}.H.calls.gz"

echo "Final calls file:"
echo "${OUT_DIR}/${REFNAME}.chr_${CHR_NUM}.H.calls.gz"

# =========================================================
# 12. CLEAN UP TEMPORARY HELPER
# =========================================================

rm -rf -- "$TMPDIR_SAMPLE"

echo "======================================"
echo "FINISHED CHROMOSOME $CHROM"
echo "======================================"
