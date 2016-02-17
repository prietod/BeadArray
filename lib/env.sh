# lib/env.sh

setvar() { var="$1"; shift; value="$@"; eval $var=\"${value}\"; }

#-------------------------------------------------------------------------
# Config
#-------------------------------------------------------------------------

setvar JINFINITI_DATA_DIR             "/hiidata/teddy/data/jinfiniti/gene_expression"
setvar RAW_DATA_DIR                   "$BASE/tmp/data/raw"
setvar FILTERED_DATA_DIR              "$BASE/tmp/data/filtered"
setvar RAW_CHIP_LIST                  "$BASE/tmp/work/common/chip-list-raw.txt"
setvar FILTERED_CHIP_LIST             "$BASE/tmp/work/common/chip-list-filtered.txt"
setvar RAW_DATA_EXCLUDE_PATTERNS      "$BASE/static/exclude-data-raw.txt"

setvar EXCLUDE_ALL_PATTERN_FILE       "$BASE/tmp/work/common/exclude-all.txt"
setvar EXCLUDE_NON_DONOR_PATTERN_FILE "$BASE/tmp/work/common/exclude-non-donor.txt"
setvar EXCLUDE_DONOR_PATTERN_FILE     "$BASE/tmp/work/common/exclude-donor.txt"

setvar METHODS                        "${METHODS:-$(echo {a..p})}"
setvar DATASETS                       "${DATASETS:-donor non-donor all}"

setvar GENERATE_PDFS                  "${GENERATE_PDFS:-no}" # yes or no

setvar DEFAULT_PIPELINE               "generate-raw-chip-list
                                       copy-raw-data
                                       raw-qc
                                       combine-qc-details
                                       sample-filter
                                       copy-filtered-data
                                       methods-step-1
                                       combine-expression
                                       generate-phenotype-details
                                       methods-step-2
                                       combine-normalized-methods"

setvar PIPELINE                       "${PIPELINE:-${DEFAULT_PIPELINE}}"

setvar STATE_DIR                      "$BASE/tmp/state"
setvar R_VERSION                      "3.2.3"
setvar TEST_MODE                      "${TEST_MODE:-false}"
setvar LOG_DIR                        "$BASE/tmp/log"
setvar WORK_DIR                       "$BASE/tmp/work"

#setvar R_LIBS                         "${HOME}/lib/R/${R_VERSION}"
setvar R_LIBS                         "${BASE}/.R/${R_VERSION}"

# Our own variables we use when running Slurm
setvar SLURM_MEM                      "16G"
setvar SLURM_STATE_DIR                "$BASE/tmp/slurm"
setvar SLURM_LOG_DIR                  "$BASE/tmp/log"

# Slurm Environmental Variables
setvar SBATCH_PARTITION               "hii02"
setvar SBATCH_TIMELIMIT               "0-08:00:00"
setvar SLURM_CPUS_PER_TASK            "2"

export \
  R_VERSION \
  R_LIBS \
  SBATCH_PARTITION \
  SBATCH_TIMELIMIT \
  SLURM_CPUS_PER_TASK
