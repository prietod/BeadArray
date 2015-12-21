#!/usr/bin/env bash

echo "MYBASE: $MYBASE"

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

echo my_process $SLURM_ARRAY_TASK_ID
