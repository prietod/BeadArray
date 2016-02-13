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
