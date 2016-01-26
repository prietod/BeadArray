#!/usr/bin/env bash

BASE="$(readlink -f $(readlink -f $(dirname "${BASH_SOURCE}")/../..))"

script=$1
expression_file=$2

find $BASE/tmp/work/${script}/chip -name "*[0-9]_${expression_file}" \
  | $BASE/bin/util/fix-data > $BASE/tmp/work/${script}/combined/${expression_file}
