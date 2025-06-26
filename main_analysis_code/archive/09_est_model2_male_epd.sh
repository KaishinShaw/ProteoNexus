#!/bin/bash
#SBATCH --job-name=estmod2
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/est_model2_male_d4_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/est_model2_male_d4_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-162%100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=48:00:00

# Set parameters
EST=/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/08_est_model2_epd.R
SEX=male
let k=0

for EL in `seq 1 54`
do
	for TL in `seq 55 57`
	do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
		then
			Rscript ${EST} --sex ${SEX} --exp_id ${EL} --out_id ${TL}
		fi
	done
done