#!/bin/bash
#SBATCH --partition=bigmem-1TB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=2-00:00
#SBATCH --job-name=simcov_test_w_model
#SBATCH --mail-user=
#SBATCH --mail-type=END

#./lungmodel 8000 8000 1 8627 7010 0
./lungmodel 25255 21031 43734 0 0 0
