#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o fst.%A_%a.out
#SBATCH --array=0-54%10

set -euo pipefail

# =========================================================
# 1. SETTINGS
# =========================================================

BASE="/scratch/leuven/357/vsc35707/popgen"

cd "$BASE"

CALLS_DIR="${BASE}/final-calls"
OUTDIR="${BASE}/final-stats"
POPSFILE="${BASE}/popfile.txt"
SCRIPT="${BASE}/popgenWindows_egglib.py"

PYTHON="/data/leuven/357/vsc35707/miniconda3/bin/python"

mkdir -p "$OUTDIR"

# =========================================================
# 2. CHROMOSOME AND POPULATION PAIR
# =========================================================

CHR_INDEX=$((SLURM_ARRAY_TASK_ID / 5 + 1))
PAIR_INDEX=$((SLURM_ARRAY_TASK_ID % 5))

pop1=(
    B_NIE_T
    F_GUE_T
    P_AVE_T
    S_HUE_T
    S_BAR_T
)

pop2=(
    B_DUD_S
    F_GUE_S
    P_AVE_S
    S_COT_S
    S_CET_S
)

POP1="${pop1[$PAIR_INDEX]}"
POP2="${pop2[$PAIR_INDEX]}"

CALLS_FILE="${CALLS_DIR}/Pchal_Bar_LW.chr_${CHR_INDEX}.H.calls.gz"

OUTFILE="${OUTDIR}/Pchal_Bar_LW.chr_${CHR_INDEX}.pop${PAIR_INDEX}.stats"

echo "======================================"
echo "Chromosome: CHR${CHR_INDEX}"
echo "Comparison: ${POP1} vs ${POP2}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "======================================"

# =========================================================
# 3. VALIDATE INPUTS
# =========================================================

for FILE in "$CALLS_FILE" "$POPSFILE" "$SCRIPT"; do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing input: $FILE" >&2
        exit 1
    fi
done

gzip -t "$CALLS_FILE"

# Confirm that the requested populations exist.
for POP in "$POP1" "$POP2"; do
    N=$(awk -v p="$POP" '$2 == p {n++} END {print n+0}' "$POPSFILE")

    if [[ "$N" -eq 0 ]]; then
        echo "ERROR: No samples assigned to $POP" >&2
        exit 1
    fi

    echo "$POP: $N individuals"
done

# =========================================================
# 4. CALCULATE WINDOW STATISTICS
# =========================================================

echo "Running population genomic statistics..."

"$PYTHON" "$SCRIPT" \
    -w 25000 \
    -s 12500 \
    --minSites 250 \
    --maxMissing 0.25 \
    -T 10 \
    --windType coordinate \
    -f phased \
    -g "$CALLS_FILE" \
    --popsFile "$POPSFILE" \
    -o "$OUTFILE" \
    -p "$POP1" \
    -p "$POP2" \
    -eggB FstWC,Dxy \
    -eggW Pi,D

# =========================================================
# 5. VALIDATE OUTPUT
# =========================================================

if [[ ! -s "$OUTFILE" ]]; then
    echo "ERROR: Output file is missing or empty." >&2
    exit 1
fi

echo "======================================"
echo "COMPLETED CHR${CHR_INDEX}"
echo "Comparison: ${POP1} vs ${POP2}"
echo "Output: $OUTFILE"
echo "======================================"
