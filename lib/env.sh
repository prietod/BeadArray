#!/usr/bin/env bash

MYDIR=$(readlink -f $(dirname $0))
MYBASE=$(readlink -f $MYDIR/..)

filter_files() {
  grep -v '\.zip$' | grep -v '/level2/' | grep -v 'map_info.txt'
}

find_cache_and_lookup_base_directory() {
  if [[ ! -f $MYBASE/tmp/filedb.txt || -n $RESET  ]]; then
    find /hiidata/teddy/data/jinfiniti/gene_expression -type f | filter_files | sort > $MYBASE/tmp/filedb.txt
  fi
}

find_vial_candidates() {
  grep '\.txt$' $MYBASE/tmp/filedb.txt | awk -F/ '{print $NF}' | sed 's/\(.*\)_[A-Z]\.txt$/\1/' | sort -u
}

find_vial_instances() {
  for c in $( find_vial_candidates ); do
    count=$( grep "/${c}_[A-Z]\.txt" $MYBASE/tmp/filedb.txt  | wc -l )
    echo $c:$count
  done
}

get_sample_txt_files() {
  find $BASEDIR -type f -name '*.txt'
}

#------------------------------------------------------------------------
# main
#------------------------------------------------------------------------

if [[ $1 == '-r' ]]; then
  RESET='y'
  shift
fi

find_cache_and_lookup_base_directory

#find_vial_candidates

find_vial_instances


#get_sample_txt_files
#bundles="201505"
#
#for bundle in $bundles; do
#  find $basedir/$bundle -type f | grep -v '\.zip'  | awk -F'/' '{print $NF}'
#done

# 3998787099_B_Grn.locs
#!/usr/bin/env bash
#

filter_excluded() {
  grep -v '/level2' | grep -v '7911755025' | grep -v '\.zip$'
}

setup_barcode_dirs() {
  local chip_barcode=$1
  work_dir=$WORK_DIR/${chip_barcode}

  work_data_dir=${work_dir}/data
  work_qc_details=${work_dir}/qc_details
  work_results_dir=${work_dir}/results

  for d in $work_data_dir $work_qc_details $work_results_dir; do
    mkdir -p $d
  done
}

#set -x

MYBASE="$( dirname $0 )/.."

DATA_DIR=${DATA_DIR:-/hiidata/teddy/data/jinfiniti/gene_expression}
WORK_DIR=${WORK_DIR:-$MYBASE/tmp}

CODE=${CODE:-$MYBASE/code/BeadArray.R}

CHIP_BARCODES="${@:-3998787078}"

module load apps/R/${R_VERSION:-3.2.3}

for chip_barcode in ${CHIP_BARCODES}; do

  setup_barcode_dirs ${chip_barcode}

  for f in $( find $DATA_DIR -type f -name "${chip_barcode}*" | filter_excluded ); do
     [[ -f ${work_data_dir}/$( basename $f ) ]] || cp -v $f ${work_data_dir}
   done

  R --slave --quiet --no-restore --no-save \
    --args ${work_data_dir} ${work_qc_details} ${work_results_dir} < ${CODE}
done

#!/usr/bin/env bash

MYDIR=$(readlink -f $(dirname $0))
MYBASE=$(readlink -f $MYDIR/..)

export MYBASE

datestamp=$( date +%FT%T )

log_dir=$MYBASE/tmp/log

(cd $log_dir && mkdir $datestamp && ln -nsf $datestamp current)



sbatch \
  --job-name=BeadArray \
  --output=$log_dir/$datestamp/%a.log \
  --array=1-2 \
  --time=00:01:00 \
  --partition=hii02 \
  --ntasks=1 \
  --wrap=$MYBASE/bin/run $1 <command string>

  $MYBASE/bin/submit.sbatch
#!/usr/bin/env bash

echo "MYBASE: $MYBASE"

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

echo my_process $SLURM_ARRAY_TASK_ID
#!/usr/bin/env bash

MYDIR=$(readlink -f $(dirname $0))
MYBASE=$(readlink -f $MYDIR/../..)

DATA_BASEDIR=${DATA_BASEDIR:-/hiidata/teddy/data/jinfiniti/gene_expression}
INDEX_FILE=${INDEX_FILE:-$MYBASE/tmp/run/filedb.txt}

filter_out_excluded_files_and_dirs() {
  egrep -v '(/level2/|\.zip$|map_info.txt$)'
}

convert_tabs_to_colons() {
  sed 's/\t/:/g'
}

strip_first_line() {
  sed -n '2,$p'
}

print_first_line() {
  head -1 | sed 's/^/#/'
}

sort_by_chip_barcode() {
  sort -t: -k7,7
}

#------------------------------------------------------------------------
# main
#------------------------------------------------------------------------

mkdir -p "$(dirname $INDEX_FILE)"

#---

echo "$(basename $0): Creating Filtered/Sorted File Index '$INDEX_FILE' for path '$DATA_BASEDIR'" 1>&2

find $DATA_BASEDIR -type f | filter_out_excluded_files_and_dirs | sort > $INDEX_FILE

#---

echo "$(basename $0): Reformatting '$DATA_BASEDIR/map_info.txt' and writing to '$MYBASE/tmp/run/map_info.csv'" 1>&2

print_first_line < $DATA_BASEDIR/map_info.txt | convert_tabs_to_colons > $MYBASE/tmp/run/map_info-header.csv

strip_first_line < $DATA_BASEDIR/map_info.txt | convert_tabs_to_colons | sort_by_chip_barcode > $MYBASE/tmp/run/map_info.csv

#!/usr/bin/env bash

MYBASE="$( dirname $0 )/../.."

if [[ $1 == '-u' ]]; then
  filter="sort -u"
  shift
else
  filter="cat"
fi

case $1 in
  subject_id) field=1;;
  vial_barcode_number) field=2;;
  Test_Name) field=3;;
  Donor_Number) field=4;;
  Box) field=5;;
  Row) field=6;;
  Chip_Barcode) field=7;;
  Array) field=8;;
  draw_dte) field=9;;
  sample_status) field=10;;
  Date_Received_Sample) field=11;;
  Date_of_Evaluation) field=12;;
  Comments) field=13;;

  *) echo "Usage: $( basename $0) <subject_id|vial_barcode_number|Test_Name|Donor_Number|Box|Row|Chip_Barcode|Array|draw_dte|sample_status|Date_Received_Sample|Date_of_Evaluation|Comments>" 1>&2; echo 1>&2; exit 1;;

esac

shift

awk -F: "{print \$$field }" $MYBASE/tmp/run/map_info.csv | $filter
