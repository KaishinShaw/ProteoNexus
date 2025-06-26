#!/bin/bash
#SBATCH --job-name=tot
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/est_tot_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/est_tot_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-549
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=1:00:00

# Set parameters
EST=/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/07_est_total_emd.R
LABEL=/gpfs/chencao/ysbioinfor/project/proteohubProject/02_data/model1_label.txt
let k=0

for ID in $(cat ${LABEL})
do
	for SEX in all female male
	do 
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
		then
			Rscript ${EST} --sex ${SEX} --fac_id ${ID}
		fi
	done
done
