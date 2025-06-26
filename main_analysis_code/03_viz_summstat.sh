#!/bin/bash
#SBATCH --job-name=viz_summ
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/viz_summ3_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/viz_summ3_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-1000%300
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=1:00:00

# Set parameters
VIZ_SUMM=/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/03_viz_summstat.R
PRO=/gpfs/chencao/ysbioinfor/Datasets/ukb/pheno/04_olink/protein_label3.txt
let k=0

for pro in $(cat ${PRO})
do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
	then
		Rscript ${VIZ_SUMM} --pro ${pro}
	fi
done
