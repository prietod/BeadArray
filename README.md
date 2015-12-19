# BeadArray

BeadArray HPC Project

## Installation

    git clone git@github.com:USF-HII/BeadArray.git 

    cd BeadArray

    bin/run

## Scratch

### Shell Scripts

All shell scripts in `bin/` should use the convention for their beginning stanza:

    #!/usr/bin/env bash

    MYDIR=$(readlink -f $(dirname $0))
    MYBASE=$(readlink -f $MYDIR/..)


### Directories

- `bin/` - Driver scripts, utility scripts, mostly if not all written in Bash.
- `code/` - The intelligent bits that do the real analysis work. 
- `tmp/` - Directory on which all temporary data is stored (e.g. data, runtime data, results, qc) - It will be ignored on Git Commits

#### Temp Directory

- `tmp/data` - Data copied into a certain directory format for analysis.
- `tmp/qc` - Quality-Control Information
- `tmp/results` - Final Results of analysis 
- `tmp/run` - Temporary scratch files, cached data, etc.

 
