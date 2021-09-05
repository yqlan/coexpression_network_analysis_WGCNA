#!/bin/sh
#############################################################################################################################
# Date: Sat Mar 13 09:48:04 CST 2021
###############################################################################################################################

cd /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/temp/rawdata/paper16_paper17
num=1
for file in `ls *fastq`
do
tag=${file%_*}
echo "
#PBS -N mp.${num}
#PBS -q core40
#PBS -l mem=20gb,walltime=99:00:00,nodes=1:ppn=2
#HSCHED -s RNAseq+hisat2+mm

ref_dir=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/reference
rawdata_dir=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/temp/rawdata/paper16_paper17
outputdir=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/mapping

mkdir \${outputdir}/$tag

export PATH=/asnas/mishl_group/lanyq/software/hisat2-2.2.1:\$PATH

# single end data mapping
hisat2 -p 2 -x \$ref_dir/GRCm39.Gchr.genome --dta -U \$rawdata_dir/$file -S \${outputdir}/$tag/$tag.sam
" > /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/code/mapping.$tag.sh
num=$[$num+1]
done


# pair end data mapping
#hisat2 -p 5 -x \$ref_dir/GRCm39.Gchr.genome --dta -1 \$rawdata_dir/example/reads/reads_1.fa -2 \$rawdata_dir/example/reads/reads_2.fa -S \${outputdir}/$tag/$tag.sam
