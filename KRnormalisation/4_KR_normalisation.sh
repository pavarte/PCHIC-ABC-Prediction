#!/bin/bash
#SBATCH --job-name=KR_normalisation         # Job name
#SBATCH --time=6:00:00                      # Time limit hrs:min:sec
#SBATCH --cpus-per-gpu=4 
#SBATCH --gres=gpu:2 
#SBATCH --mem=400G
#SBATCH --output=KR_normalisation%j.out
#SBATCH --error=KR_normalisation%j.err
#SBATCH --partition=gpu


module load miniconda3
conda activate octave
cd ./

for i in $(ls subsampled);
do
        echo $i
        gcMapExplorer normKR -i ${i} -fi ccmap -o ${i}.KR.normalised.ccmap -fo ccmap -m RAM
        python3 KRnorm_vector.py --ccmap_input ${i} --output_file_path .
        python3 ccamp_to_text_export.py --ccmap_input ${i}.KR.normalised.ccmap --output ${i}.KRobserved.txt
done
