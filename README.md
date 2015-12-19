# BeadArray

BeadArray HPC Project

## Installation

    git clone git@github.com:USF-HII/BeadArray.git 

    cd BeadArray

    bin/run

## Scratch

### Shell Scripts

All shell scripts in `bin/` should use the following convention at the top of their file: 

    #!/usr/bin/env bash

    MYDIR=$(readlink -f $(dirname $0))
    MYBASE=$(readlink -f $MYDIR/..)


### Directories

- `bin/` - Driver scripts, utility scripts, mostly if not all written in Bash.
- `code/` - The intelligent bits that do the real analysis work. 
- `tmp/` - Directory on which all temporal data is stored (input data, utility data, qc, results) - It will be ignored upon any Git Commit and should be copied out as needed. 

#### tmp/ 

- `tmp/data` - Data copied into a certain directory format for analysis.
- `tmp/qc` - Quality-Control Information
- `tmp/results` - Final Results of analysis 
- `tmp/run` - Temporary scratch files, cached data, etc.

 
