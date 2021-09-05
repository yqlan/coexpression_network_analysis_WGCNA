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

#PBS -N ev.st.${num}
#PBS -q core40
#PBS -l mem=10gb,walltime=99:00:00,nodes=1:ppn=2
#HSCHED -s RNAseq+stringtie+mm

genes_gtf=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/reference/gencode.vM26.annotation.Gchr.gtf

bamdir=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/mapping

outputdir=/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/expression_value

export PATH=/asnas/mishl_group/lanyq/software/stringtie-2.1.5:\$PATH

mkdir \${outputdir}/$samp

stringtie -p 2 -G \$genes_gtf -B -A \${outputdir}/$samp/$samp.gene_abund.txt -e -o \${outputdir}/$samp/$samp.out.gtf \$bamdir/$samp/$samp.bam

" > /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/code/fpkm_stringtie.$samp.sh

num=$[$num+1]

done

