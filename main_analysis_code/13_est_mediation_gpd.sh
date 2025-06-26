#!/bin/bash
#SBATCH --job-name=med
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/est_mediation_gpd_all_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/est_mediation_gpd_all_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-401%150
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=24:00:00

# Set parameters
EST=/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/13_est_mediation_gpd.R
SEX=all
let k=0

for SN in `seq 1 401`
do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
	then
		Rscript ${EST} --sex ${SEX} --snp_num ${SN}
	fi
done