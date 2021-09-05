#!/bin/sh
#############################################################################################################################
# Date: Thu Mar 11 16:02:17 CST 2021
#
###############################################################################################################################
#PBS -N sra2fq
#PBS -q core40
#PBS -l mem=20gb,walltime=99:00:00,nodes=1:ppn=5
#HSCHED -s RNAseq+fastqdump+mm

cd /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/temp/rawdata/paper16_paper17


export LD_LIBRARY_PATH=/software/biosoft/software/sra-tools-master/download/tools/lib64:$LD_LIBRARY_PATH
export NCBI_VDB_LIBDIR=/software/biosoft/software/sra-tools-master/download/tools/lib64
export PATH=/software/biosoft/software/sra-tools-master/biosoft/bin:$PATH

for file in `ls *.sra`
do
fastq-dump ${file}
done

