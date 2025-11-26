#!/bin/bash
# --- Configuration ---
PROJECT_ROOT_DIR=$(pwd)
LOG_SUBDIR="logs"
R_SCRIPT_RELATIVE_PATH="scripts/S02b_fit_alfak_vary_misseg_sims.R"   # Path to your R script from project root
INPUT_FILES_RELATIVE_DIR="data/raw/ABM_varyMisseg"                    # Path to input files from project root

LOG_DIR="${PROJECT_ROOT_DIR}/${LOG_SUBDIR}"
R_SCRIPT_FULL_PATH="${PROJECT_ROOT_DIR}/${R_SCRIPT_RELATIVE_PATH}"
INPUT_FILES_FULL_DIR="${PROJECT_ROOT_DIR}/${INPUT_FILES_RELATIVE_DIR}"

mkdir -p "${LOG_DIR}"

INPUT_FILES_BASENAME_ARRAY=()

if [ ! -d "${INPUT_FILES_FULL_DIR}" ] || [ ! -r "${INPUT_FILES_FULL_DIR}" ]; then
    echo "Error: Input directory ${INPUT_FILES_FULL_DIR} not found or not readable. Exiting."
    exit 1
fi

while IFS= read -r -d $'\0' file_path; do
    INPUT_FILES_BASENAME_ARRAY+=("$(basename "$file_path")")
done < <(find "${INPUT_FILES_FULL_DIR}" -maxdepth 1 -name "*.Rds" -print0)

if [ ${#INPUT_FILES_BASENAME_ARRAY[@]} -eq 0 ]; then
    echo "Error: No .Rds files found in ${INPUT_FILES_FULL_DIR}. Exiting."
    exit 1
fi

NUM_TOTAL_TASKS="${#INPUT_FILES_BASENAME_ARRAY[@]}"

SLURM_JOB_NAME_PREFIX=${0##*/} && SLURM_JOB_NAME_PREFIX=${SLURM_JOB_NAME_PREFIX%%.*}

echo "Submitting ${NUM_TOTAL_TASKS} tasks for ${SLURM_JOB_NAME_PREFIX} to SLURM."
echo "R_SCRIPT: ${R_SCRIPT_FULL_PATH}"
echo "LOGS_DIR: ${LOG_DIR}"

sbatch <<EOSUBMIT
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-06:00:00
#SBATCH --qos=small
#SBATCH --array=0-$((${NUM_TOTAL_TASKS}-1))
#SBATCH --job-name=${SLURM_JOB_NAME_PREFIX}
#SBATCH --output=${LOG_DIR}/${SLURM_JOB_NAME_PREFIX}_%A_%a.out
#SBATCH --error=${LOG_DIR}/${SLURM_JOB_NAME_PREFIX}_%A_%a.err

# --- SLURM Task ---
# Ensure your R environment is set up (e.g., module load R/version) if Rscript is not in default PATH
ml R/4.2.0

TASK_FILES_LIST_EXPANDED=(${INPUT_FILES_BASENAME_ARRAY[*]})

CURRENT_INPUT_FILE_BASENAME="\${TASK_FILES_LIST_EXPANDED[\$SLURM_ARRAY_TASK_ID]}"
INPUT_FILES_PARENT_DIR="${INPUT_FILES_FULL_DIR}"
FULL_PATH_TO_CURRENT_INPUT_FILE="\${INPUT_FILES_PARENT_DIR}/\${CURRENT_INPUT_FILE_BASENAME}"

echo "SLURM Array Task ID: \${SLURM_ARRAY_TASK_ID} | File: \${CURRENT_INPUT_FILE_BASENAME}"

# Usage: Rscript name_of_this_script.R <input_file>
R_COMMAND="Rscript ${R_SCRIPT_FULL_PATH} \"\${FULL_PATH_TO_CURRENT_INPUT_FILE}\""
eval "\${R_COMMAND}"
EXIT_STATUS=\$?

if [ \$EXIT_STATUS -ne 0 ]; then
    echo "Error in SLURM Array Task ID: \${SLURM_ARRAY_TASK_ID} | File: \${CURRENT_INPUT_FILE_BASENAME} | Exit Status: \${EXIT_STATUS}"
fi
exit \${EXIT_STATUS}
EOSUBMIT

echo "SLURM job array submitted. Check status with squeue -u \$USER and logs in ${LOG_DIR}"


