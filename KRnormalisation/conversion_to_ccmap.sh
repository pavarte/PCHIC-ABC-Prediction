#!/bin/bash
#SBATCH --job-name=KR_normalisation         # Job name
#SBATCH --time=6:00:00                      # Time limit hrs:min:sec
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --mem=50G                          # Memory limit
#SBATCH --output=KR_normalisation%j.out
#SBATCH --partition=nice
#SBATCH --ntasks=1

module load miniconda3
conda activate gcmap_env
cd ./

gcMapExplorer coo2cmap -i input_raw_list.txt -ccm RawObserved_normal -r 5kbp -od /gc_expl_test
