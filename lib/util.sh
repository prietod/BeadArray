# general functions

sbatch-submit() {
  local job_name=$1; shift
  local slurm_opts=$1; shift
  local script=$1; shift

  if echo "${slurm_opts}" | grep -q -- '--array'; then
    ensure_dir ${SLURM_LOG_DIR}/${job_name}
    local log_spec="${SLURM_LOG_DIR}/${job_name}/%a.log"
  else
    local log_spec="${SLURM_LOG_DIR}/${job_name}.log"
  fi

  [[ -f ${script} ]] || error_exit "Script '${script}' not found."
  [[ -x ${script} ]] || error_exit "Script '${script}' not executable."

  job_id=$( sbatch \
             --job-name=${job_name} \
             --output=${log_spec} \
             --mail-type=FAIL \
             ${slurm_opts} \
             ${script} "$@" | awk '{print $NF}' )

  echo ${job_id}
}

run() {
  local name=$1; shift

  ensure_dir ${STATE_DIR}

  if ! [[ -f ${STATE_DIR}/${name} ]]; then
    echo
    echo "-----------------------------------------------------------------------------"
    echo "run> ${name} "$@""
    echo "-----------------------------------------------------------------------------"
    echo

    eval "$name "$@""

    [[ $? -eq 0 ]] && touch ${STATE_DIR}/${name}
  fi
}

error_exit() {
  local msg=$1
  echo "Error: $msg" | tee $FAILED_DIR/${script_short}/${chip_barcode}.txt 1>&2
  exit 1
}

mark_job_complete() {
  date +%FT%T > $COMPLETE_DIR/${script_short}/${chip_barcode}.txt
}

check_have_group() {
  local group=$1
  local id_output=$( id )
  if ! echo ${id_output} | grep -q "\<${group}\>"; then
    error "Exiting. Node dropped '${group}' group ownership in session, data inaccessible. id-output: ${id_output}"
  fi
}

error() {
  echo "$@" 1>&2
}

error_exit() {
  error "Error $(basename $0): "$@""
  exit 1
}

ensure_dir() { local path=$1; [[ -d ${path} ]] || mkdir -p ${path}; }

line_count() { wc -l "$1" | awk '{print $1}'; }

join_lines() {
  tr "\n" ":" | sed 's/:$//'
}

get_job_ids() {
  local job_name
  for job_name in "$@"; do
    cat ${SLURM_STATE_DIR}/${job_name}-job_id.txt
  done | join_lines
}
