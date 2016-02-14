# BeadArray HPC Pipeline

## Installation

    git clone git@github.com:USF-HII/BeadArray.git

## Run the pipeline

    $ ./bin/run-all

## Google R/Shell Style Guides

- https://google.github.io/styleguide/Rguide.xml
- https://google.github.io/styleguide/shell.xml

## Directory Structure

- `bin/` - Driver scripts, utility scripts, mostly if not all written in Bash.

- `code/` - The intelligent bits that do the real analysis work.

- `tmp/work/<job_name>/` - Work directory for each job.

- `tmp/log/<job_name>.log` - Slurm log for each job.

- `tmp/log/<job_name>/<n>.log` - Slurm logs for each task in an array job.

## Notes

combine-expression()
  sbatch/combine-expression
    echo ./tmp/work/method-a-step-1-all/chip/9250939009/results/9250939009_control_expression_lumi.txt \
      | ./bin/util/combine-expression-data lumi

- install-R-libs"
- generate-raw-chip-list
- copy-raw-data
- run-raw-qc
- combine-qc-details
- sample-filter
- copy-filtered-data
- methods-step-1
- combine-expression > sbatch/combine-expression-submit -> sbatch/combine-expression
- phenotype-details
- methods-step-2

Work-around:

    ran method-{a..p}-step-1-donor and its complete while non-donor and all finish step-1

    $ ls -l ./tmp/work/method-a-step-1-donor/chip/*/results

    $ export DATASETS="donor"

    $ bin/run-all combine-expression

    $ bin/util/find-chip-arrays tmp/data/filtered/donor/ > tmp/work/common/chip-list-filtered-donor.txt

    $ wc -l tmp/work/common/chip-list-filtered-donor.txt
    66 tmp/work/common/chip-list-filtered-donor.txt

    $ head -1 tmp/work/common/BeadArray_phenotype_details-all.txt \
        > tmp/work/common/BeadArray_phenotype_details-donor.txt

    $ grep -f tmp/work/common/chip-list-filtered-donor.txt tmp/work/common/BeadArray_phenotype_details-all.txt \
        >> tmp/work/common/BeadArray_phenotype_details-donor.txt


