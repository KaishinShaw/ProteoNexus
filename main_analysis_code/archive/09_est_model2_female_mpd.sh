#!/bin/bash
#SBATCH --job-name=estmod2
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/est_model2_female_mpd_d1_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/est_model2_female_mpd_d1_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-969%200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=48:00:00

# Set parameters
EST=/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/09_est_model2_mpd.R
SEX=female
let k=0

for ML in `seq 1 17`
do
	for TL in `seq 1 57`
	do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
		then
			echo ${ML}
			echo ${TL}
			Rscript ${EST} --sex ${SEX} --ms_id ${ML} --out_id ${TL}
		fi
	done
done
