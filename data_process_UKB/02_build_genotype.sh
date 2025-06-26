#!/bin/bash

PLINK=/public/home/biostat03/biosoft/plink1.9
EID=/public/home/Datasets/ukb/geno/EUR_protien/eid

chr=22

INFILE=/public/home/Datasets/ukb/geno/EUR/hm3/chr${chr}
OUTFILE=/public/home/Datasets/ukb/geno/EUR_protien/chr${chr}

${PLINK} --bfile ${INFILE} --keep ${EID} --indiv-sort f ${EID}
