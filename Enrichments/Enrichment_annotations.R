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
compute_enrichments=function(snps0.01,snps0.001,snps0.0001){
  
  ## enrichments
  Numbers=NULL
  coeffs=NULL
  for(i in 1:length(annotations)){
    annot=annotations[[i]]
    annot_sel0.01=nrow((data %>% filter((Conseq %in% annot)  & (SNP %in% snps0.01) & transcript_type=="protein_coding")))
    annot_NotSel=nrow((data %>% filter((Conseq %in% annot) & !(SNP %in% snps0.01) & transcript_type=="protein_coding")))
    NotAnnot_sel=nrow((data %>% filter(!(Conseq %in% annot)   & (SNP %in% snps0.01) & transcript_type=="protein_coding")))
    NotAnnot_NotSel=nrow((data %>% filter(!(Conseq %in% annot) & !(SNP %in% snps0.01) & transcript_type=="protein_coding")))
    res0.01=oddsratio(matrix(c(annot_sel0.01,annot_NotSel,NotAnnot_sel,NotAnnot_NotSel),byrow=T,nrow=2))
    
    
    annot_sel0.001=nrow((data %>% filter((Conseq %in% annot) & (SNP %in% snps0.001) & transcript_type=="protein_coding")))
    annot_NotSel=nrow((data %>% filter((Conseq %in% annot) & !(SNP %in% snps0.001) & transcript_type=="protein_coding")))
    NotAnnot_sel=nrow((data %>% filter(!(Conseq %in% annot) & (SNP %in% snps0.001) & transcript_type=="protein_coding")))
    NotAnnot_NotSel=nrow((data %>% filter(!(Conseq %in% annot) & !(SNP %in% snps0.001) & transcript_type=="protein_coding")))
    res0.001=oddsratio(matrix(c(annot_sel0.001,annot_NotSel,NotAnnot_sel,NotAnnot_NotSel),byrow=T,nrow=2))
    
    annot_sel0.0001=nrow((data %>% filter((Conseq %in% annot) & (SNP %in% snps0.0001) & transcript_type=="protein_coding")))
    annot_NotSel=nrow((data %>% filter((Conseq %in% annot) & !(SNP %in% snps0.0001) & transcript_type=="protein_coding")))
    NotAnnot_sel=nrow((data %>% filter(!(Conseq %in% annot) & (SNP %in% snps0.0001) & transcript_type=="protein_coding")))
    NotAnnot_NotSel=nrow((data %>% filter(!(Conseq %in% annot) & !(SNP %in% snps0.0001) & transcript_type=="protein_coding")))
    res0.0001=oddsratio(matrix(c(annot_sel0.0001,annot_NotSel,NotAnnot_sel,NotAnnot_NotSel),byrow=T,nrow=2))
    
    coeffs=rbind(coeffs,c(as.character(as.vector(names[i])),
                          res0.01$measure[2,1],res0.01$measure[2,2],res0.01$measure[2,3],
                          res0.001$measure[2,1],res0.001$measure[2,2],res0.001$measure[2,3],
                          res0.0001$measure[2,1],res0.0001$measure[2,2],res0.0001$measure[2,3]))
    Numbers=rbind(Numbers,c(annot_sel0.01,annot_sel0.001,annot_sel0.0001))
    cat(i,"\n")
  }
  Altogether=cbind(coeffs,Numbers)
  Altogether=data.frame(Altogether)
  colnames(Altogether)=c("Annotation","OR_1%","OR_0.025_1%","OR_0.975_1%","OR_0.1%","OR_0.025_0.1%","OR_0.975_0.1%",
                         "OR_0.01%","OR_0.025_0.01%","OR_0.975_0.01%","N_1%","N_0.1%","N_0.01%")
  for(i in 2:ncol(Altogether)){
    Altogether[,i]=as.numeric(as.vector(Altogether[,i]))
  }
  
  ## formatting
  boxLabels = as.character(as.vector(Altogether[,1]))
  df <- data.frame(labels=rep(boxLabels,3),
                   boxOdds = c(Altogether[,2],Altogether[,5],Altogether[,8]), 
                   boxCILow = c(Altogether[,3],Altogether[,6],Altogether[,9]), 
                   boxCIHigh = c(Altogether[,4],Altogether[,7],Altogether[,10]),
                   N = c(Altogether[,11],Altogether[,12],Altogether[,13]),
                   model = c(rep(paste("p < 0.01",sep=""),length(boxLabels)),
                             rep(paste("p < 0.001",sep=""),length(boxLabels)),
                             rep(paste("p < 1e-04",sep=""),length(boxLabels))))
  
  df$boxOdds=as.numeric(as.vector(df$boxOdds))
  df$boxCILow=as.numeric(as.vector(df$boxCILow))
  df$boxCIHigh=as.numeric(as.vector(df$boxCIHigh))
  df$N=as.numeric(as.vector(df$N))
  df$model=as.factor(as.vector(df$model))
  df$labels=as.character(as.vector(df$labels))
  
  boxLabels2=rep(boxLabels,4)
  df$model <- factor(df$model, levels = c(paste("p < 0.01",sep=""),
                                          paste("p < 0.001",sep=""),
                                          paste("p < 1e-04",sep="")))
  df$labels <- factor(df$labels, levels = boxLabels)
  df
}


### FIGURE 2A

### ANNOTATION ENRICHMENTS

## They were conducted on variants annotated in protein_coding transcripts (removes intergenic variants and allows
## a more homogeneous view of the impact of a given annotation)

datum=(data %>% filter(pemp<1e-02 & (GroupsLD %in% groups) & transcript_type=="protein_coding"))

## possible consequences on translation
conseq=(data %>% filter(transcript_type=="protein_coding") )$Conseq

# annotations to study
missense=names(table(conseq))[grepl("^missense",names(table(conseq)))]
UTR=names(table(conseq))[grepl("prime_UTR",names(table(conseq)))]
Upstream_Downstream=names(table(conseq))[grepl("downstream_gene_variant|upstream_gene_variant",names(table(conseq)))]
intron="intron_variant"
synonymous="synonymous_variant"
annotations=list(missense,UTR,Upstream_Downstream,intron,synonymous)
names=c("Missense","UTR","Upstream/Downstream","Intron","Synonymous")


## function to compute enrichments
## first within candidate LD groups
snps0.01=(data %>% filter(pemp<1e-02 & (GroupsLD %in% groups)))$SNP
snps0.001=(data %>% filter(pemp<1e-03 & (GroupsLD %in% groups)))$SNP
snps0.0001=(data %>% filter(pemp<1e-04 & (GroupsLD %in% groups)))$SNP
within_groupsLD=compute_enrichments(snps0.01,snps0.001,snps0.0001)

## secondly, across all loci
snps0.01=(data %>% filter(pemp<1e-02))$SNP
snps0.001=(data %>% filter(pemp<1e-03))$SNP
snps0.0001=(data %>% filter(pemp<1e-04))$SNP
across_loci=compute_enrichments(snps0.01,snps0.001,snps0.0001)

## put together to plot
df=rbind(within_groupsLD,across_loci)
df$Data=c(rep("89 loci",nrow(within_groupsLD)),rep("All loci",nrow(across_loci)))
df$Data=factor(df$Data,levels=c("All loci","89 loci"))
df$labels=factor(df$labels,levels=rev(c("Missense","Synonymous","UTR","Upstream/Downstream","Intron")))
df=data.frame(df)

## RESULT
df

