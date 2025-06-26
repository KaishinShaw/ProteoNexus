#!/bin/bash
#SBATCH --job-name=med
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/im_est_mediation_emd_female_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/im_est_mediation_emd_female_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-108
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=24:00:00

# Set parameters
EST=/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/09_est_mediation_emd.R
SEX=female
GROUP=3
let k=0

for PP in `seq 1 108`
do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
	then
		echo ${ML}
		echo ${TL}
		Rscript ${EST} --sex ${SEX} --pair ${PP} --group ${GROUP}
	fi
done