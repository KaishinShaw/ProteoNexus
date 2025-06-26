#!/bin/bash
#SBATCH --job-name=mk_summ
#SBATCH --error=/gpfs/chencao/ysbioinfor/project/proteohubProject/err_file/mk_summ3_%a.err
#SBATCH --out=/gpfs/chencao/ysbioinfor/project/proteohubProject/out_file/mk_summ3_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-1000%100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cu,privority,batch01,gpu,fat
#SBATCH --time=12:00:00

# Set parameters
GENO_PATH=/gpfs/chencao/ysbioinfor/Datasets/ukb/geno/EUR_protein/hm3/
PRO=/gpfs/chencao/ysbioinfor/Datasets/ukb/pheno/04_olink/protein_label3.txt
GEMMA=/gpfs/chencao/ysbioinfor/software/gemma-0.98.1-linux-static
let k=0

for pro in $(cat ${PRO})
do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
	then
		for sex in male female
		# for sex in all
		do
			cd /gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl/${pro}/
			INFILE=${GENO_PATH}${sex}/merge
			${GEMMA} -bfile ${INFILE} -n ${k} -notsnp -lm 1 -o summ_${sex}
			rm -rf ./output/summ_${sex}.log.txt
			gzip -f ./output/summ_${sex}.assoc.txt
		done
	fi
done
