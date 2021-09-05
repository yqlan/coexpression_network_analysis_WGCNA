#!/bin/sh
#############################################################################################################################
# Date: Fri Mar 12 20:28:34 CST 2021
###############################################################################################################################
#PBS -N index-build
#PBS -q q512G
#PBS -l mem=200gb,walltime=99:00:00,nodes=1:ppn=2
#HSCHED -s RNAseq+hisat2+mm

cd /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/reference
export PATH=/asnas/mishl_group/lanyq/software/hisat2-2.2.1:$PATH
hisat2-build -p 2 --ss spliceSites.txt --exon exons.txt GRCm39.Gchr.genome.fa GRCm39.Gchr.genome
