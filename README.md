# bEADaRRay

## Installation

    git clone git@github.com:USF-HII/BeadArray.git

## Running Locally

    bin/run <R_script> [chip_barcode]

For example:

    bin/run code/BeadArray_qc.R 3998755068

## Running via Slurm

To run for all chip_barcodes:

   bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/Bead_Array_qc.R

## Google R Style Guide

- https://google.github.io/styleguide/Rguide.xml

## R Code Arguments

The `bin/run` script will copy the necessary files into a working data structure and
will pass 2 arguments to your R code.

Your R code should contain the following lines:

```
#------------------------------------------------------------------------
# Get cli args and assign appropriate variables
#------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

raw_data_dir <- toString(args[1])
raw_qc_dir <- toString(args[2])
```

At this point you can reference these variables to locate the desired directory for reading/writing.

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

Use the command `bin/util/get-mapinfo` to query this file and use the flag `-u` to list a column uniquely.

## Directory Structure

- `bin/` - Driver scripts, utility scripts, mostly if not all written in Bash.

- `code/` - The intelligent bits that do the real analysis work.

- `tmp/` - Directory in which all temporal data is stored (input data, utility data, qc, results).
           It will be ignored upon any Git Commit and should be copied out as needed.

### tmp/

    - tmp/
      - data/
         - <chip_barcode>/...
           - <data_files>...
      - work
         - <code_name>
           - combined/
              - <combined_files>
           - chips/
             - <chip-barcode>/...
                - raw/
                  - qc/
                    - <generated_files...>

## Notes

### Chip Bar Code 9224522100

We are experiencing a segmentation fault on chip barcode `9224522100` on our HPC platform
but were able to analyze this on a different system.

After the `Beadarray_qc` and `Beadarray_qc_average` pipelines are run we must unarchive the results
for this chip as follows:

From the project base directory of the project:

    for code_name in BeadArray_qc BeadArray_qc_average; do
      dir=tmp/work/${code_name}/chips/9224522100/raw/qc
      [[ -d ${dir} ]] && rm -rf ${dir}
      mkdir -p ${dir}
      tar --directory ${dir} -zxf /hiidata/projects/BeadArray/${code_name}-chips-9224522100-raw-qc.tar.gz
    done

## Scratch

### High Level Flow

- Apply `static/exclude-data-raw.txt` and copy all `/hiidata/teddy/data/jinfiniti/gene_expression` list to `tmp/data/raw`

- For `BeadArray_{qc,qc_average}.R` as `<script>`

    - Multiple run `code/<script>` on `tmp/data/raw/<chip>` generating `tmp/work/<script>/chips/<chip>/raw/qc/<files..>`

- Concatenate all `tmp/work/<script>/chips/<chip>/raw/qc/<chip>_details.txt` to `tmp/work/<script>/combined/raw-qc-details.txt`

- Single run `code/BeadArray_sample_filter.R` on `tmp/work/<script>/combined/raw-qc-details.txt`

- Verify both `tmp/work/BeadArray_{qc,qc_average}.R/combined/raw-qc-details.txt` and copy first
  column to `tmp/work/all/exclude-data-filtered.txt`

- Apply `tmp/work/all/exclude-data-fitlered.txt` and copy data from `tmp/data/raw/<chip>` to `tmp/data/filtered`

- Multipe run `code/BeadArray_method_1_step_{1..5}.R` on `tmp/data/filtered/<chip>`

- Multiple run `code/BeadArray_approach_a_step_1.R` on `tmp/work/<script>/chips/<chip>/raw-qc-details.txt`


### Thoughts

#### Left-Oriented versus Right-Oriented Namespace Variance

General idea is if you vary the left-hand side of a namespace (e.g. foo/bar/...),
you should not vary the right-hand side of the directories or files (e.g. foo/bar/foo-bar-abc.txt)
as it this double-ended variance is not conducive to automated pipelines.

[need example showing why its not conducive ;-)]

The reason this double-ended variance exists is because one forsees in the future that
the two namespaces will somehow be combined and a conflict will occur.

With automated pipelines, it should be trivial to to apply the Left-Oriented variance, E.g.:

    foo/bar/abc.txt
    foo/baz/abc.txt

when consolidation occurs:

    combined/foo-bar-abc.txt
    combined/foo-baz-abc.txt

