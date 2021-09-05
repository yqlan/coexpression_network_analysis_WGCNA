#!/bin/sh
#############################################################################################################################
# Date: Sat Mar 13 09:48:04 CST 2021
###############################################################################################################################

cd /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/temp/rawdata
num=17
for file in `ls *1.fastq.gz`
do
tag=${file%_*}
echo "
#PBS -N mp.${num}
#PBS -q core40
#PBS -l mem=20gb,walltime=99:00:00,nodes=1:ppn=5
#HSCHED -s RNAseq+hisat2+mm

ref_dir=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/reference
rawdata_dir=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/temp/rawdata
outputdir=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/mapping

mkdir \${outputdir}/$tag

export PATH=/asnas/mishl_group/lanyq/software/hisat2-2.2.1:\$PATH

# pair end data mapping
hisat2 -p 5 -x \$ref_dir/GRCm39.Gchr.genome --dta -1 \$rawdata_dir/${tag}_1.fastq.gz -2 \$rawdata_dir/${tag}_2.fastq.gz -S \${outputdir}/$tag/$tag.sam

" > /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/code/mapping.$tag.sh
num=$[$num+1]
done


