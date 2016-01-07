# BeadArray

BeadArray HPC Project

## Installation

    git clone git@github.com:USF-HII/BeadArray.git

    cd BeadArray

    bin/run

## Info

### Shell Scripts

All shell scripts in `bin/` should use the following convention at the top of their file:

    #!/usr/bin/env bash

    MYDIR=$(readlink -f $(dirname $0))
    MYBASE=$(readlink -f $MYDIR/..)


### Directories

- `bin/` - Driver scripts, utility scripts, mostly if not all written in Bash.
- `code/` - The intelligent bits that do the real analysis work.
- `tmp/` - Directory in which all temporal data is stored (input data, utility data, qc, results) - It will be ignored upon any Git Commit and should be copied out as needed.

#### tmp/

- `tmp/data` - Data extracted and/or copied into structure for analysis.
- `tmp/qc` - Quality-Control Information
- `tmp/results` - Final Results of analysis
- `tmp/run` - Temporary scratch files, cached data, etc.


## Scratch

tmp/
  123456/
    data/ # (hiidata)
    qc_details/
      *.pdf
      *.txt
      qc/
        <f1>/
        <f2>/
        <f3>/...
     results/
       ...


   post-processing-script
    - roll-up data
    - generate aggregate image



tmp/1234456/data <---
tmp/1234456/qc <--- pdf
tmp/1234456/results

----

    for n in $( bin/get-mapinfo Chip_Barcode -u ); do
      grep ":$n:" tmp/run/map_info.csv | sort -t: -k8,8 ;
      echo
    done > tmp/run/map_info-by-chip_barcode.txt

    for n in $( bin/get-mapinfo subject_id -u ); do
      grep "^$n:" tmp/run/map_info.csv | sort -t: -k8,8 ;
      echo
    done > tmp/run/map_info-by-subject_id.txt


