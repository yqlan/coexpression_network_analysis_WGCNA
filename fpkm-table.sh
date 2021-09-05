#/bin/sh

cd /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/expression_value

for samp in `ls -d *`; do awk -v OFS="\t" '(NR>1){if(!fpkm[$1";"$2]){fpkm[$1";"$2]=$8}else{fpkm[$1";"$2]=fpkm[$1";"$2]+$8}}END{for (k in fpkm){print k,fpkm[k]}}' ./$samp/$samp.gene_abund.txt | sort -k 1,1 | awk -v OFS="\t" -v samp=$samp 'BEGIN{print "geneid_genename","fpkm_"samp}{print $0}' > ./$samp/$samp.FPKM.txt; done

paste */*.FPKM.txt | awk '{printf $1;for(i=1;i<=NF;i++){if(i%2==0){printf "\t"$i}};printf "\n"}' > fpkm-table.txt

# genes expressed in 25% samples (7s), and retain only coding and lncRNA (including TEC) genes. Data then is log() transformed.
awk -v OFS="\t" '(NR>1){m=0;for(i=2;i<=NF;i=i+1){if($i >= 0.1){m++}};if(m>6){for(i=2;i<=NF;i=i+1){$i=log($i+1)};print $0}}' fpkm-table.txt | sort -k 1,1 | join -t $'\t' - /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/reference/gencode.vM26.lncRNA-codingene.id.name.txt | awk -v OFS="\t" 'NR==1{print $0}NR>FNR{print $0}' fpkm-table.txt - > logfpkm_repeat25sample.txt

# genes expressed in 30% samples (8s), and retain only coding and lncRNA (including TEC) genes. Data then is log() transformed.
awk -v OFS="\t" '(NR>1){m=0;for(i=2;i<=NF;i=i+1){if($i >= 0.1){m++}};if(m>7){for(i=2;i<=NF;i=i+1){$i=log($i+1)};print $0}}' fpkm-table.txt | sort -k 1,1 | join -t $'\t' - /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/reference/gencode.vM26.lncRNA-codingene.id.name.txt | awk -v OFS="\t" 'NR==1{print $0}NR>FNR{print $0}' fpkm-table.txt - > logfpkm_repeat25sample_8s.txt

# genes expressed in 25% samples (6s), and retain only coding and lncRNA (including TEC) genes. Data then is log() transformed.
awk -v OFS="\t" '(NR>1){m=0;for(i=2;i<=NF;i=i+1){if($i >= 0.1){m++}};if(m>5){for(i=2;i<=NF;i=i+1){$i=log($i+1)};print $0}}' fpkm-table.txt | sort -k 1,1 | join -t $'\t' - /asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/reference/gencode.vM26.lncRNA-codingene.id.name.txt | awk -v OFS="\t" 'NR==1{print $0}NR>FNR{print $0}' fpkm-table.txt - > logfpkm_repeat25sample_6s.txt
