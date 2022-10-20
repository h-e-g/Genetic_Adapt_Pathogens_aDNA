library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)

### Figure 4
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
data_neg=datas[[datatypes[2]]]
colnames(data)[2]="CHR.POS"
colnames(data_neg)[2]="CHR.POS"
colnames(data)[3]="CHR.POS.COUNTED.ALT"
colnames(data_neg)[3]="CHR.POS.COUNTED.ALT"

data=data %>% filter(!(abs(Epoch5.F-Epoch7.F)>0.1))
data_neg=data_neg %>% filter(!(abs(Epoch7.F-Epoch5.F)>0.1))

## Read immunity genes
iig=as.character(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GO_ENRICHMENT/IIG.txt",header=F)[,1]))

## LD groups under positive selection to remove from the analysis of negative selection (including groups + classic loci under pos sel)
group_to_remove_neg=loadRData("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ESTIMATION_WORKSPACE/tmp/group_to_remove_neg")

### FUNCTIONS
compute_enrichments=function(snps0.01,snps0.001,snps0.0001){
  
  ## enrichments
  Numbers=NULL
  coeffs=NULL
  for(i in 1:length(annotations)){
    annot=annotations[[i]]
    annot_sel0.01=nrow((data_neg %>% filter((Conseq %in% annot)  & (SNP %in% snps0.01) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    annot_NotSel=nrow((data_neg %>% filter((Conseq %in% annot) & !(SNP %in% snps0.01) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    NotAnnot_sel=nrow((data_neg %>% filter(!(Conseq %in% annot)   & (SNP %in% snps0.01) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    NotAnnot_NotSel=nrow((data_neg %>% filter(!(Conseq %in% annot) & !(SNP %in% snps0.01) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    res0.01=oddsratio(matrix(c(annot_sel0.01,annot_NotSel,NotAnnot_sel,NotAnnot_NotSel),byrow=T,nrow=2))
    
    
    annot_sel0.001=nrow((data_neg %>% filter((Conseq %in% annot) & (SNP %in% snps0.001) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    annot_NotSel=nrow((data_neg %>% filter((Conseq %in% annot) & !(SNP %in% snps0.001) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    NotAnnot_sel=nrow((data_neg %>% filter(!(Conseq %in% annot) & (SNP %in% snps0.001) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    NotAnnot_NotSel=nrow((data_neg %>% filter(!(Conseq %in% annot) & !(SNP %in% snps0.001) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    res0.001=oddsratio(matrix(c(annot_sel0.001,annot_NotSel,NotAnnot_sel,NotAnnot_NotSel),byrow=T,nrow=2))
    
    annot_sel0.0001=nrow((data_neg %>% filter((Conseq %in% annot) & (SNP %in% snps0.0001) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    annot_NotSel=nrow((data_neg %>% filter((Conseq %in% annot) & !(SNP %in% snps0.0001) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    NotAnnot_sel=nrow((data_neg %>% filter(!(Conseq %in% annot) & (SNP %in% snps0.0001) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
    NotAnnot_NotSel=nrow((data_neg %>% filter(!(Conseq %in% annot) & !(SNP %in% snps0.0001) & transcript_type=="protein_coding" & !(GroupsLD %in% group_to_remove_neg))))
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

## TABLE S6
##### enrichment of missense variants among negatively selected variants
conseq=(data_neg %>% filter(transcript_type=="protein_coding") )$Conseq

# annotations to study
missense=names(table(conseq))[grepl("^missense",names(table(conseq)))]
UTR=names(table(conseq))[grepl("prime_UTR",names(table(conseq)))]
Upstream_Downstream=names(table(conseq))[grepl("downstream_gene_variant|upstream_gene_variant",names(table(conseq)))]
intron="intron_variant"
synonymous="synonymous_variant"
annotations=list(missense,UTR,Upstream_Downstream,intron,synonymous)
names=c("Missense","UTR","Upstream/Downstream","Intron","Synonymous")

snps0.01=(data_neg %>% filter(pemp<1e-02))$SNP
snps0.001=(data_neg %>% filter(pemp<1e-03))$SNP
snps0.0001=(data_neg %>% filter(pemp<1e-04))$SNP
outside_positive=compute_enrichments(snps0.01,snps0.001,snps0.0001)



##### enrichment of low-frequent variants
## Test the enrichment of variants below a given percentile (1%, 0.1%)
## P values were defined for bins of frequency so for each bin of frequency, we only expect 1% or 0.1% of variants below the 1% or 0.1% percentile, respectively.

sig=c(1e-02,1e-04)
thres=c(seq(0,0.05,0.025),0.1,0.2,0.8)
Enrichment=NULL
for(j in 2:(length(thres)-1)){
  for(i in 1:length(sig)){
    NbrSigBin1=nrow((data_neg %>% filter(!(GroupsLD %in% group_to_remove_neg) & MAX>thres[j] & MAX<=thres[j+1] & pemp<sig[i])))
    NbrBin1=nrow((data_neg %>% filter(!(GroupsLD %in% group_to_remove_neg) & MAX>thres[j] & MAX<=thres[j+1] & pemp>=sig[i])))
    NbrSigExpBin1=round(nrow((data_neg %>% filter(!(GroupsLD %in% group_to_remove_neg) & MAX>thres[j] & MAX<=thres[j+1])))*(sig[i]))
    NbrExpBin1=nrow((data_neg %>% filter(!(GroupsLD %in% group_to_remove_neg) & MAX>thres[j] & MAX<=thres[j+1])))-NbrSigExpBin1
    res=oddsratio(matrix(c(NbrSigBin1,NbrSigExpBin1,NbrBin1,NbrExpBin1),nrow=2,byrow=T))
    Enrichment=rbind(Enrichment,c(format(res$measure[2,1],digits=2),format(res$measure[2,2],digits=2),format(res$measure[2,3],digits=2),res$p.value[2,3],paste(thres[j],thres[j+1],sep="-"),paste("Model p <",sig[i],sep=""),NbrSigBin1,"Negative"))
    
  }
}
colnames(Enrichment)=c("OR","OR_l","OR_u","P","Bin","model","N","Type")
