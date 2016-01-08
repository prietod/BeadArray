# BeadArray

## Installation

    git clone git@github.com:USF-HII/BeadArray.git
    cd BeadArray

## Running Locally

    bin/run <R_script> [chip_barcode] [chip_barcode...]

For example:

    bin/run code/BeadArray.R 3998755068

## Running via Slurm

Generate all the bar codes:

   bin/util/get-mapinfo -u chip_barcode > tmp/chip_barcodes.txt

   bin/slurm-submit code/BeadArray.R < tmp/chip_barcodes.txt

## Env File

The file `lib/env.sh` contains environmental variables which have defaults but which may be overriden.

At the time of this documentation:

```
export DATA_DIR=${DATA_DIR:-/hiidata/teddy/data/jinfiniti/gene_expression}
export MAP_INFO=${MAP_INFO:-/hiidata/teddy/data/jinfiniti/gene_expression/map_info.txt}
```

## Map Info

The map info file is a tab-delimited list with the following fields:

- subject_id
- vial_barcode_number
- test_name
- donor_number
- box
- row
- chip_barcode
- array
- draw_dte
- sample_status
- date_received_sample
- date_of_evaluation
- comments

Use the command `bin/util/get-mapinfo` to query this file.

## Style Guide

- https://google.github.io/styleguide/Rguide.xml

## R Code Arguments

The `bin/run` script will copy the necessary files into a working data structure and will pass 3 arguments to your R code.

Your R code should contain the following lines:

```
#------------------------------------------------------------------------
# Get cli arguments
#------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

data_dir <- toString(args[1])
results_dir <- toString(args[2])
qc_details_dir <- toString(args[3])
```

At this point you can reference these variables to locate the desired directory for reading/writing.

## Working Directory structure

The default Working Directory `$WORK_DIR` is under the project at `tmp/work` and contains the following sub-directories
for each `chip_barcode`:

- `data_dir/` - contains all of the files for processing for a chip_barcode
- `qc_details_dir/` - all necessary QC details :-)
- `results_dir/` - all results

For example:

    - tmp/work/3998755068/data/<files..>
    - tmp/work/3998755068/results/<files..>
    - tmp/work/3998755068/qc_details/<files..>

Example of `data_dir` contents for one vial_barcode:

    tmp/work/3998755068/data/3998755068_A.csv
    tmp/work/3998755068/data/3998755068_A.txt
    tmp/work/3998755068/data/3998755068_A_Grn.idat
    tmp/work/3998755068/data/3998755068_A_Grn.locs
    tmp/work/3998755068/data/3998755068_A_Grn.tif
    tmp/work/3998755068/data/3998755068_A_Grn.xml
    tmp/work/3998755068/data/3998755068_B.csv ...etc.

### Directories

- `bin/` - Driver scripts, utility scripts, mostly if not all written in Bash.
- `code/` - The intelligent bits that do the real analysis work.
- `tmp/` - Directory in which all temporal data is stored (input data, utility data, qc, results) - It will be ignored upon any Git Commit and should be copied out as needed.


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



        for n in $( bin/get-mapinfo Chip_Barcode -u ); do
          grep ":$n:" tmp/run/map_info.csv | sort -t: -k8,8 ;
          echo
        done > tmp/run/map_info-by-chip_barcode.txt

        for n in $( bin/get-mapinfo subject_id -u ); do
          grep "^$n:" tmp/run/map_info.csv | sort -t: -k8,8 ;
          echo
        done > tmp/run/map_info-by-subject_id.txt


