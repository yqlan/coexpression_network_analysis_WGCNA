#!/bin/sh
#############################################################################################################################
# Date: Sat Mar 13 09:48:04 CST 2021
###############################################################################################################################
cd /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/mapping
num=1
for samp in `ls -d *`
do

echo "
#!/bin/sh

#PBS -N sort.${num}
#PBS -q core40
#PBS -l mem=20gb,walltime=99:00:00,nodes=1:ppn=2
#HSCHED -s RNAseq+samtools+mm

export PATH=/software/biosoft/software/samtools1.9/bin:\$PATH

cd /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/mapping

samtools view -Su ./$samp/$samp.sam | samtools sort - -T ./$samp/$samp -o ./$samp/$samp.bam

" > /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/code/sortSam2Bam.samtools.$samp.sh

num=$[$num+1]

done
