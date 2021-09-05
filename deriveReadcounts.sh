#! /bin/sh

cd /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/expression_value

for file in `ls | grep -v txt`;

do 

rl=`samtools view ../mapping/$file/$file.bam | head -1 | awk '{print length($10)}'`;

python /asnas/mishl_group/lanyq/software/stringtie-2.1.5/prepDE.py -g ./$file/$file-gene_count_matrix.csv -t ./$file/$file-transcript_count_matrix.csv -l $rl -p $file

done


python /asnas/mishl_group/lanyq/software/stringtie-2.1.5/prepDE.py -g ../diffExprAnaly/csh3_no4HT-4HT/edgeR/csh3.gene_count_matrix.csv -t ../diffExprAnaly/csh3_no4HT-4HT/edgeR/csh3.transcript_count_matrix.csv -l 101 -p csh3_

python /asnas/mishl_group/lanyq/software/stringtie-2.1.5/prepDE.py -g ../diffExprAnaly/nascentRNAseq_0-24and0-48/edgeR/nascentRNAseq.gene_count_matrix.csv -t ../diffExprAnaly/nascentRNAseq_0-24and0-48/edgeR/nascentRNAseq.transcript_count_matrix.csv -l 49 -p nascentRNAseq_

python /asnas/mishl_group/lanyq/software/stringtie-2.1.5/prepDE.py -i ciGdal1-scr.txt -g ../diffExprAnaly/rnaseq_CiGdal/edgeR/ciGdal-scr.gene_count_matrix.csv -t ../diffExprAnaly/rnaseq_CiGdal/edgeR/ciGdal-scr.transcript_count_matrix.csv -l 150

