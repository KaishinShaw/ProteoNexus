PROJ_PATH=/gpfs/chencao/ysbioinfor/project/proteohubProject
PLINK=/gpfs/chencao/ysbioinfor/software/plink1.9

for SEX in male all 
do
	INFILE=/gpfs/chencao/ysbioinfor/Datasets/ukb/geno/EUR_protein/hm3/${SEX}/merge
	SIG_SNP=${PROJ_PATH}/02_data/sig_genotype/sig_snp_${SEX}.txt
	OUTFILE=${PROJ_PATH}/02_data/sig_genotype/sig_${SEX}
	${PLINK} --bfile ${INFILE} --extract ${SIG_SNP} --recode A --out ${OUTFILE}
done