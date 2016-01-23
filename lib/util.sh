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
