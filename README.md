# BeadArray

BeadArray HPC Project

## Installation

    git clone git@github.com:USF-HII/BeadArray.git 

    cd BeadArray

    bin/run

## Scratch

### Shell Scripts

All shell scripts in `bin/` should use the convention:

  #!/usr/bin/env bash

  MYDIR=$(readlink -f $(dirname $0))
  MYBASE=$(readlink -f $MYDIR/..)

for their beginning stanza.

### Directories

- `tmp/` - Directory on which all temporary data is stored (e.g. data, runtime data, results, qc) - It will be ignored on Git Commits

 
