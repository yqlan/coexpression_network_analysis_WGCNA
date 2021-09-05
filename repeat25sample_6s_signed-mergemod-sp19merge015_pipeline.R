

getwd();
setwd("/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/coexpression");
library(WGCNA);
options(stringsAsFactors = FALSE);
enableWGCNAThreads();
###########################################   < 1 >    #################################################################################
files <- dir();
for (expfile in files[grep("logfpkm_repeat25sample_lnc-pc.txt", files)]){	## names of the first colum must be geneid_genename  ##
genexpr<-read.table(expfile,header=TRUE,sep="\t");
dim(genexpr);
genexpr_0 = as.data.frame(t(genexpr[,-1]));
dim(genexpr_0);
names(genexpr_0)<-genexpr$geneid_genename;
rownames(genexpr_0)<-names(genexpr)[-1];
dim(genexpr_0);
genexpr <- genexpr_0;
dim(genexpr);
fnl<-strsplit(expfile,"_");
save(genexpr,file = paste(paste(fnl[[1]][1],fnl[[1]][2],sep="_"),"RData",sep="."));
}

###########################################   < 2 >    #################################################################################
files <- dir();
for (rdata in files[grep("logfpkm_repeat25sample.RData", files)]){                  ######   
data<-load(rdata);
sampleTree <- hclust(dist(genexpr), method = "average");
fn<-sub(pattern="\\.RData", replacement="",rdata);
pdf(file = paste(fn,"SampleClustering.pdf",sep="_"), width = 12, height = 9);
par(cex = 1);
par(mar = c(2,5,3,2));
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()
}

powers = c(5:25);
files <- dir();
for (rdata in files[grep("logfpkm_repeat25sample.RData", files)]){                ######
data<-load(rdata);
sft <- pickSoftThreshold(genexpr, powerVector = powers, networkType = "signed", verbose = 5);
fn<- sub(pattern="logfpkm_",replacement="",sub(pattern="\\.RData", replacement="",rdata));
pdf(file = paste(fn,"soft-powers.pdf",sep="_"), width = 9, height = 5);
par(mfrow = c(1,2));
cex1 <- 0.7;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.5,col="red");
abline(h=0.8,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"));
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");
dev.off();
}

##########################################     < 2.TOM >    ############################################################################
data<-load("logfpkm_repeat25sample.RData");          #####
softPower <- 26;    ### R.sq>=0.80 ###
adjacency <- adjacency(genexpr, type = "signed", power = softPower);    ##singed/unsigned##
save(adjacency,file = "repeat25sample_adjacency.RData");    ## save adjacency ##                      ######
TOM <- TOMsimilarity(adjacency,TOMType = "signed");    ##singed/unsigned##
dissTOM <- 1-TOM;
geneTree <- hclust(as.dist(dissTOM), method = "average");
minModuleSize <- 30;
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize);
table(dynamicMods);
dynamicColors <- labels2colors(dynamicMods);
table(dynamicColors);
pdf(file = "repeat25sample_clustersOnTOM-basedDissimilarityAndGene-modules.pdf", width=8, height=6);          ######
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Modules",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "clusters on TOM-based dissimilarity and gene-modules");
dev.off();
save(TOM,file = "repeat25sample_TOM.RData");                   ######
save(softPower,geneTree,dynamicMods,dynamicColors,file = "repeat25sample_sftpower_dynamicLables_dynamicColors_geneTree.RData");          #####


MEList <- moduleEigengenes(genexpr, colors = dynamicColors, softPower = softPower,);
MEs <- MEList$eigengenes;
MEDiss <- 1-cor(MEs);
METree <- hclust(as.dist(MEDiss), method = "average");
pdf(file = "repeat25sample_ClusteringOfModuleEigengenes.pdf",width=7, height=6);           ######
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "");
MEDissThres = 0.15;	####
abline(h=MEDissThres, col = "red");
dev.off();


merge = mergeCloseModules(genexpr,dynamicColors,cutHeight = MEDissThres,verbose = 3);
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
pdf(file = "repeat25sample_clustersOnTOM-basedDissimilarityAndGene-modules_mergedModules.pdf", wi = 12, he = 9);         ######
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged Modules"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05);
dev.off();
moduleColors = mergedColors;
colorOrder = c("grey", standardColors(100));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
save(MEs, moduleLabels, moduleColors, file = "repeat25sample_MEs_ModuleLables_ModuleColors_geneTree.RData");           ######

##################################################    < 3 >    ########################################################################
data1<-load("logfpkm_repeat25sample.RData");  #####
data2<-load("repeat25sample_sftpower_dynamicLables_dynamicColors_geneTree.RData");  ####
data3<-load("repeat25sample_MEs_ModuleLables_ModuleColors_geneTree.RData")    ####
datTraits <- read.table("repeat25sample_leukstate.txt",header=TRUE);   ####

nGenes <- ncol(genexpr);
nSamples <- nrow(genexpr);
MEs0 <- moduleEigengenes(genexpr, moduleColors, softPower = softPower)$eigengenes;
MEs <- orderMEs(MEs0);
moduleTraitCor <- cor(MEs, datTraits);
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples);
pdf(file="repeat25sample_Module-trait_relationships.pdf", width = 8, height = 15);  ####
textMatrix <- paste(signif(moduleTraitCor,2),"(",signif(moduleTraitPvalue,1),")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(3, 8, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,cex.lab=0.6,xLabels = names(datTraits),yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,colors = greenWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.4,zlim = c(-1,1),main = paste("Module-trait relationships"));
dev.off();

modNames <- substring(names(MEs),3);
geneModuleMembership = as.data.frame(cor(genexpr, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", names(MEs), sep="");
names(MMPvalue) = paste("p.MM", names(MEs), sep="");
geneTraitSignificance = as.data.frame(cor(genexpr, datTraits));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(datTraits), sep="");
names(GSPvalue) = paste("p.GS.", names(datTraits), sep="");

geneInfo0 <- data.frame(geneid_genename = names(genexpr),moduleColor = moduleColors,geneTraitSignificance,GSPvalue);
dim(geneInfo0);
modOrder = order(-abs(cor(MEs, datTraits)));
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0);
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),paste("p.MM.", modNames[modOrder[mod]], sep=""));
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.leuk_state));
geneInfo = geneInfo0[geneOrder, ];
write.table(geneInfo,file="repeat25sample_geneInfo.txt",row.names=F,col.names=T,sep="\t",quote=FALSE);   ####


###########################################   < 4 >    #################################################################################
annot <- read.table(file = "genelist.txt",header=TRUE);		#### genelist
probes <- names(genexpr);
probes2annot <- match(probes, annot$geneid_genename);
allgenename <- annot$genename[probes2annot];
intModules <- c("darkmagenta");   #### module name
for (module in intModules)
{
modprobes <- (moduleColors==module);
modgenes <- allgenename[modprobes];
fileName <- paste("repeat25sample_modgenes-", module, ".txt", sep="");   ####
write.table(as.data.frame(modgenes), file = fileName,row.names = FALSE, col.names = FALSE,quote=FALSE)
}

############################################   < 5 >   #################################################################################
MEs <- moduleEigengenes(genexpr, moduleColors, softPower = softPower)$eigengenes;
leuksta <- as.data.frame(datTraits$leuk_state);
names(leuksta) <- "leuksta";
MET <- orderMEs(cbind(MEs, leuksta));
pdf(file = "repeat25sample_EigengeneNetworks.pdf", width = 10, height = 16);    ####
par(cex = 0.5);
plotEigengeneNetworks(MET, "", marDendro = c(0,4,3,4.6), marHeatmap = c(4,5,0,2), cex.lab = 0.8, xLabelsAngle= 90);
dev.off();


############################################   < 6 >   #################################################################################
getwd();
setwd("/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/coexpression");		####
library("WGCNA");
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

data4<-load("repeat25sample_TOM.RData");    ####
data5<-load("logfpkm_repeat25sample.RData");    ####
data6<-load("repeat25sample_MEs_ModuleLables_ModuleColors_geneTree.RData");    ####
modules<-c();    #### module name
annot <- read.table(file = "genelist.txt",header=TRUE);
probes<-names(genexpr);
for (mod in modules) {
        inModule <- (moduleColors==mod);
        modProbes <- probes[inModule];
        modTOM <- TOM[inModule, inModule];
        dimnames(modTOM) <- list(modProbes,modProbes);
        nTop <- 30;
        IMConn <- softConnectivity(genexpr[, modProbes],type = "signed",power= softPower);    ####
        top <- (rank(-IMConn) <= nTop);
        vis_top30 = exportNetworkToVisANT(modTOM[top, top],file = paste("repeat25sample_VisANTInput-", mod, "-top30.txt", sep=""),weighted = TRUE,threshold = 0,probeToGene = data.frame(annot$geneid_genename, annot$genename));            ####
        vis = exportNetworkToVisANT(modTOM,file = paste("repeat25sample_VisANTInput-", mod, ".txt", sep=""),weighted = TRUE,threshold = 0,probeToGene = data.frame(annot$geneid_genename, annot$genename))};                 ####


data4<-load("repeat25sample_TOM.RData");    ####
data5<-load("logfpkm_repeat25sample.RData");    ####  genexpr  #####
data6<-load("repeat25sample_MEs_ModuleLables_ModuleColors_geneTree.RData");    ####  MEs, moduleLabels, moduleColors   ####
annot <- read.table(file = "genelist.txt",header=TRUE);
modules<-c("brown4");    ####  module name
probes<-names(genexpr);
inModule <- is.finite(match(moduleColors, modules));
modProbes <- probes[inModule];
modGenes <- annot$genename[match(modProbes, annot$geneid_genename)];
modTOM <- TOM[inModule, inModule];
dimnames(modTOM) <- list(modProbes,modProbes);
#  nodeAttr <- read.table("lncF73module_annotation.txt",header=F) 
#  names(nodeAttr)<-c("geneid","genetype","moduleColor","degree")
cyt <- exportNetworkToCytoscape(modTOM,edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),weighted = TRUE,threshold = 0.02,nodeNames = modProbes,altNodeNames = modGenes,nodeAttr = nodeAttr);



############################################   < 7 >   #################################################################################
getwd();
setwd("/asnas/mishl_group/lanyq/work/leukemia/RNAseq/mllenl_rnaseq_result/rnaseqRepeat/coexpression");		####
library("WGCNA");
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

data7<-load("repeat25sample_adjacency.RData")   ## load adjacency ##
data8<-load("repeat25sample_MEs_ModuleLables_ModuleColors_geneTree.RData");    ## load mergemod ##

intramoduleConn <- intramodularConnectivity(adjMat=adjacency, colors=moduleColors, scaleByMax = FALSE);
save(intramoduleConn,file = "repeat25sample_mergemod_intramoduleConn.RData");           ####
intramdconn<-data.frame(geneid_genename=rownames(intramoduleConn),intramoduleConn[,1:4]);
write.table(intramdconn,file="repeat25sample_mergemod_intramoduleConn.txt",row.names=F,col.names=T,sep="\t",quote=FALSE);           ####

intramoduleConn <- intramodularConnectivity(adjMat=adjacency, colors=moduleColors, scaleByMax = TRUE);
save(intramoduleConn,file = "repeat25sample_mergemod_intramoduleConn_scaleByMax.RData");          ####
intramdconn<-data.frame(geneid_genename=rownames(intramoduleConn),intramoduleConn[,1:4]);
write.table(intramdconn,file="repeat25sample_mergemod_intramoduleConn_scaleByMax.txt",row.names=F,col.names=T,sep="\t",quote=FALSE);         ####

