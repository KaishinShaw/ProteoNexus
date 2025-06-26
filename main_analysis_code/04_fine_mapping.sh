#!/bin/bash
#SBATCH --job-name=fine_map
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/fine_map_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/fine_map_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-549
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=1:00:00

# Set parameters
FINE_MAP=/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/04_fine_mapping.R
let k=0

for pro_g in `seq 1 973`
do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
	then
		Rscript ${FINE_MAP} --pro_g ${pro_g}
	fi
done
