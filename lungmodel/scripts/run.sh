#!/bin/bash
set -e
rm slurm*.out
rm airway.csv
make
sbatch simcov.sh