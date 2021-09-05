
########################################################  dynamic  ####################################################################################

setwd("/leofs/mishl_group/lanyq/temp/coexpnetwork/enlwt16sample_signed`");    ## wkdir ##
ginfo <- read.table(file = "repeat25sample_geneInfo_dynamicMod.txt",header=TRUE);    ## geneinfo ##
modules<-c(names(table(ginfo$dynamicColor)));
for (mod in modules) {
	inModule <- (ginfo$dynamicColor==mod);
	p_mm<-paste("p.MM.",mod,sep="");
	mm<-paste("MM.",mod,sep="");
	module_gene <- data.frame(geneid_genename =ginfo$geneid_genename[inModule],dynamicColor=ginfo$dynamicColor[inModule],GS.leuk_state=ginfo$GS.leuk_state[inModule],p.GS.leuk_state=ginfo$p.GS.leuk_state[inModule],ginfo[[mm]][inModule],ginfo[[p_mm]][inModule]);
	colnames(module_gene)[5] <- mm;
	colnames(module_gene)[6] <- p_mm;
	write.table(module_gene,file=paste("geneModule_gs-mm_", mod, ".txt", sep=""),row.names=F,col.names=T,sep="\t",quote=FALSE);}

write.table(table(ginfo$dynamicColor),file="repeat25sample_dynamicmodules.txt",row.names=F,col.names=T,sep="\t",quote=FALSE);

paste repeat25sample_dynamic_intramoduleConn.txt repeat25sample_dynamic_intramoduleConn_scaleByMax.txt | awk -v OFS="\t" '{if($1==$6){print $1,$2,$3,$8,$4,$5}}' > repeat25sample_connectivity.txt

cat geneModule_gs-mm_* | sort -k 1,1 |awk '{if(NR<16757){print $0}}' > repeat25sample_gs-mm.txt   ## gene number ##

sort -k 1,1 repeat25sample_connectivity.txt | join - repeat25sample_gs-mm.txt | awk -v OFS="\t" '{print $1,$3,$4,$7,$8,$10}' > repeat25sample_intramoduleconne-gs-mm.txt

join repeat25sample_intramoduleconne-gs-mm.txt ../../lncRNA.id > repeat25sample_intramoduleconne-gs-mm_lncRNA.txt

join repeat25sample_intramoduleconne-gs-mm.txt ../../cancergene_ong-tsg-cosmic_leukgene-inhouse_leukgenedb.id > repeat25sample_intramoduleconne-gs-mm_ong-tsg-cosmic_leukgene-inhouse_leukgenedb.txt

join repeat25sample_intramoduleconne-gs-mm.txt ../../mllenl-DEgene.id > repeat25sample_intramoduleconne-gs-mm_mllenl-DEgene.txt

join repeat25sample_intramoduleconne-gs-mm.txt ../../leukegene.id > repeat25sample_intramoduleconne-gs-mm_leukgene.txt


###################################################  merge  #########################################################################################

ginfo <- read.table(file = "repeat25sample_geneInfo.txt",header=TRUE);    ## geneinfo ##
modules<-c(names(table(ginfo$moduleColor)));
for (mod in modules) {
	inModule <- (ginfo$moduleColor==mod);
	p_mm<-paste("p.MM.",mod,sep="");
	mm<-paste("MM.",mod,sep="");
	module_gene <- data.frame(geneid_genename =ginfo$geneid_genename[inModule],moduleColor=ginfo$moduleColor[inModule],GS.leuk_state=ginfo$GS.leuk_state[inModule],p.GS.leuk_state=ginfo$p.GS.leuk_state[inModule],ginfo[[mm]][inModule],ginfo[[p_mm]][inModule]);
	colnames(module_gene)[5] <- mm;
	colnames(module_gene)[6] <- p_mm;
	write.table(module_gene,file=paste("mergeModule_gs-mm_", mod, ".txt", sep=""),row.names=F,col.names=T,sep="\t",quote=FALSE);}

write.table(table(ginfo$moduleColor),file="repeat25sample_mergemodules.txt",row.names=F,col.names=T,sep="\t",quote=FALSE);

paste repeat25sample_mergemod_intramoduleConn.txt repeat25sample_mergemod_intramoduleConn_scaleByMax.txt | awk -v OFS="\t" '{if($1==$6){print $1,$2,$3,$8,$4,$5}}' > repeat25sample_mergemod_connectivity.txt

cat mergeModule_gs-mm* | sort -k 1,1 |awk '{if(NR<17062){print $0}}' > repeat25sample_mergemod_gs-mm.txt   ## gene number is 17061 ##

sort -k 1,1 repeat25sample_mergemod_connectivity.txt | join - repeat25sample_mergemod_gs-mm.txt | awk -v OFS="\t" '{print $1,$3,$4,$7,$8,$10}' > repeat25sample_mergemod_intramoduleconne-gs-mm.txt

join repeat25sample_mergemod_intramoduleconne-gs-mm.txt lncRNA.id > repeat25sample_mergemod_intramoduleconne-gs-mm_lncRNA.txt	# lncRNA id #

join -v 1 repeat25sample_mergemod_intramoduleconne-gs-mm.txt lncRNA.id > repeat25sample_mergemod_intramoduleconne-gs-mm_codingene.txt

grep darkorange repeat25sample_mergemod_intramoduleconne-gs-mm.txt | awk '{split($1,a,";");split(a[1],b,".");print b[1],$0}' > repeat25sample_mergemod_intramoduleconne-gs-mm_darkorange.txt	### module name of Gdal1: darkorange ###

join -1 2 repeat25sample_mergemod_intramoduleconne-gs-mm_darkorange.txt lncRNA.id | awk '{print $1,"lncRNA","darkorange",$4}' > lncF73module_lncRNA.txt

join -1 2 -v 1 repeat25sample_mergemod_intramoduleconne-gs-mm_darkorange.txt lncRNA.id | awk '{print $1,"mRNA","darkorange",$4}' > lncF73module_mRNA.txt

cat lncF73module_mRNA.txt lncF73module_lncRNA.txt | awk -v OFS="\t" '{split($1,a,";");split(a[1],b,".");print b[1],$2,$3,$4}' |sort -k 1,1 > lncF73module_annotation.txt

#join repeat25sample_mergemod_intramoduleconne-gs-mm.txt ../../cancergene_ong-tsg-cosmic_leukgene-inhouse_leukgenedb.id > repeat25sample_mergemod_intramoduleconne-gs-mm_ong-tsg-cosmic_leukgene-inhouse_leukgenedb.txt

#join repeat25sample_mergemod_intramoduleconne-gs-mm.txt ../../mllenl-DEgene.id > repeat25sample_mergemod_intramoduleconne-gs-mm_mllenl-DEgene.txt

#join repeat25sample_mergemod_intramoduleconne-gs-mm.txt ../../leukegene.id > repeat25sample_mergemod_intramoduleconne-gs-mm_leukgene.txt
