#!/bin/bash
#SBATCH --job-name=med
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/est_mediation_emd_male_miss_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/est_mediation_emd_male_miss_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-131
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=36:00:00

# Set parameters
EST=/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/check/09_est_mediation_emd.R
SEX=male
let k=0

for PP in `seq 1 131`
do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
	then
		echo ${ML}
		echo ${TL}
		Rscript ${EST} --sex ${SEX} --pair ${PP} 
	fi
done