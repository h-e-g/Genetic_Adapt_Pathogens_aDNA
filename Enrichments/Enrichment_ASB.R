library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggcorrplot)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
datatypes=c("DERIVED_STEP10","DERIVED_STEP30","POSITIVE_DERIVED_STEP10","POSITIVE_DERIVED_STEP30")

### LOADING DATA

FILE_NAME="REAL_RESULTS_CAPTURE_LIST_PLUS_NOT_ANNOT_NOINDEL_LDGROUPS"
datas=loadRData(paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ESTIMATION_WORKSPACE/tmp/",FILE_NAME,sep=""))
data=datas[[datatypes[4]]]
data=data %>% filter(!(abs(Epoch5.F-Epoch7.F)>0.1))
groups=loadRData("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ESTIMATION_WORKSPACE/tmp/groupsLDEnrichment0.075pval0.01")

### FUNCTIONS
enrichment_ASB=function(snps0.01,snps0.001,snps0.0001){
  Numbers=NULL
  coeffs=NULL
  annot_sel0.01=nrow((data %>% filter((SNP %in% ASB_TF) & (SNP %in% snps0.01) )))
  annot_NotSel=nrow((data %>% filter((SNP %in% ASB_TF) & !(SNP %in% snps0.01) )))
  NotAnnot_sel=nrow((data %>% filter(!(SNP %in% ASB_TF) & (SNP %in% snps0.01) )))
  NotAnnot_NotSel=nrow((data %>% filter(!(SNP %in% ASB_TF) & !(SNP %in% snps0.01) )))
  res0.01=oddsratio(matrix(c(annot_sel0.01,annot_NotSel,NotAnnot_sel,NotAnnot_NotSel),byrow=T,nrow=2))
  
  snps0.001=(data %>% filter(pemp<1e-03  ))$SNP
  annot_sel0.001=nrow((data %>% filter((SNP %in% ASB_TF) & (SNP %in% snps0.001) )))
  annot_NotSel=nrow((data %>% filter((SNP %in% ASB_TF) & !(SNP %in% snps0.001) )))
  NotAnnot_sel=nrow((data %>% filter(!(SNP %in% ASB_TF) & (SNP %in% snps0.001) )))
  NotAnnot_NotSel=nrow((data %>% filter(!(SNP %in% ASB_TF) & !(SNP %in% snps0.001) )))
  res0.001=oddsratio(matrix(c(annot_sel0.001,annot_NotSel,NotAnnot_sel,NotAnnot_NotSel),byrow=T,nrow=2))
  
  snps0.0001=(data %>% filter(pemp<1e-04 ))$SNP
  annot_sel0.0001=nrow((data %>% filter((SNP %in% ASB_TF) & (SNP %in% snps0.0001) )))
  annot_NotSel=nrow((data %>% filter((SNP %in% ASB_TF) & !(SNP %in% snps0.0001) )))
  NotAnnot_sel=nrow((data %>% filter(!(SNP %in% ASB_TF) & (SNP %in% snps0.0001) )))
  NotAnnot_NotSel=nrow((data %>% filter(!(SNP %in% ASB_TF) & !(SNP %in% snps0.0001) )))
  res0.0001=oddsratio(matrix(c(annot_sel0.0001,annot_NotSel,NotAnnot_sel,NotAnnot_NotSel),byrow=T,nrow=2))
  
  coeffs=rbind(coeffs,c(as.character(as.vector("ALL_TF")),
                        res0.01$measure[2,1],res0.01$measure[2,2],res0.01$measure[2,3],
                        res0.001$measure[2,1],res0.001$measure[2,2],res0.001$measure[2,3],
                        res0.0001$measure[2,1],res0.0001$measure[2,2],res0.0001$measure[2,3]))
  Numbers=rbind(Numbers,c(annot_sel0.01,annot_sel0.001,annot_sel0.0001))
  
  
  Altogether=cbind(coeffs,Numbers)
  Altogether=data.frame(Altogether)
  colnames(Altogether)=c("Annotation","OR_1%","OR_0.025_1%","OR_0.975_1%","OR_0.1%","OR_0.025_0.1%","OR_0.975_0.1%",
                         "OR_0.01%","OR_0.025_0.01%","OR_0.975_0.01%","N_1%","N_0.1%","N_0.01%")
  for(i in 2:ncol(Altogether)){
    Altogether[,i]=as.numeric(as.vector(Altogether[,i]))
  }
  
  ### Get desired output to plot
  
  boxLabels = c("ASB","ASB","ASB")
  
  
  df <- data.frame(labels=rep(boxLabels,1),
                   boxOdds = c(Altogether[,2],Altogether[,5],Altogether[,8]), 
                   boxCILow = c(Altogether[,3],Altogether[,6],Altogether[,9]), 
                   boxCIHigh = c(Altogether[,4],Altogether[,7],Altogether[,10]),
                   N = c(Altogether[,11],Altogether[,12],Altogether[,13]),
                   model = c(rep(paste("p < 0.01",sep=""),1),
                             rep(paste("p < 0.001",sep=""),1),
                             rep(paste("p < 1e-04",sep=""),1)))
  
  df$boxOdds=as.numeric(as.vector(df$boxOdds))
  df$boxCILow=as.numeric(as.vector(df$boxCILow))
  df$boxCIHigh=as.numeric(as.vector(df$boxCIHigh))
  df$N=as.numeric(as.vector(df$N))
  df$model=as.factor(as.vector(df$model))
  df$labels=as.character(as.vector(df$labels))
  
  df$model <- factor(df$model, levels = c(paste("p < 0.01",sep=""),
                                          paste("p < 0.001",sep=""),
                                          paste("p < 1e-04",sep="")))
  
  df
}

### FIGURE 2A
## LOAD ASBs
ASB_TF=as.character(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ADASTRA/release_Zanthar/release_dump/All_ASB_FDR0.01_TF_rs_unique.txt")[,1]))

## function to compute enrichments
## first within candidate LD groups
snps0.01=(data %>% filter(pemp<1e-02 & (GroupsLD %in% groups)))$SNP
snps0.001=(data %>% filter(pemp<1e-03 & (GroupsLD %in% groups)))$SNP
snps0.0001=(data %>% filter(pemp<1e-04 & (GroupsLD %in% groups)))$SNP
within_groupsLD=enrichment_ASB(snps0.01,snps0.001,snps0.0001)

## secondly, across all loci
snps0.01=(data %>% filter(pemp<1e-02))$SNP
snps0.001=(data %>% filter(pemp<1e-03))$SNP
snps0.0001=(data %>% filter(pemp<1e-04))$SNP
across_loci=enrichment_ASB(snps0.01,snps0.001,snps0.0001)

## put together to plot
df=rbind(within_groupsLD,across_loci)
df$Data=c(rep("89 loci",nrow(within_groupsLD)),rep("All loci",nrow(across_loci)))
df$Data=factor(df$Data,levels=c("All loci","89 loci"))
df=data.frame(df)

## RESULT
df



