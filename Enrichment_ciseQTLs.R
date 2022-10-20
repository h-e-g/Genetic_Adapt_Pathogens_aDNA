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
enrichment_eqtl=function(snps0.01,snps0.001,snps0.0001){
  coeffs_groups=NULL
  Numbers=NULL
  for(BonferroniLevel in c(5e-8,5e-15,5e-50)){
    if(BonferroniLevel=="5e-08"){
      annot="eQTL_8"
    }else if(BonferroniLevel=="5e-15"){
      annot="eQTL_15"
    }else{
      annot="eQTL_50"
    }
    
    ## get eQTLs with the given threshold
    my_eqtls=unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter(SNP %in% (eqtls %>% filter(BonferroniP<BonferroniLevel))$SNP))$SNP)
    
    eqtls_0.01=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter((SNP %in% my_eqtls) & (SNP %in% snps0.01)))$SNP))
    non_eqtls_0.01=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter(!(SNP %in% my_eqtls) & (SNP %in% snps0.01)))$SNP))
    eqtls_non_0.01=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter((SNP %in% my_eqtls) & !(SNP %in% snps0.01)))$SNP))
    non_eqtls_non_0.01=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter(!(SNP %in% my_eqtls) & !(SNP %in% snps0.01)))$SNP))
    matrix=matrix(c(non_eqtls_non_0.01,non_eqtls_0.01,eqtls_non_0.01,eqtls_0.01),nrow=2,byrow=T)
    rows=c("non-eqtls","eqtls")
    cols=c("p<0.01","p>=0.01")
    dimnames(matrix) <- list("eQTL" = rows, "Selection" = cols)
    res0.01=oddsratio(matrix)
    
    eqtls_0.001=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter((SNP %in% my_eqtls) & (SNP %in% snps0.001)))$SNP))
    non_eqtls_0.001=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter(!(SNP %in% my_eqtls) & (SNP %in% snps0.001)))$SNP))
    eqtls_non_0.001=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter((SNP %in% my_eqtls) & !(SNP %in% snps0.001)))$SNP))
    non_eqtls_non_0.001=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter(!(SNP %in% my_eqtls) & !(SNP %in% snps0.001)))$SNP))
    matrix=matrix(c(non_eqtls_non_0.001,non_eqtls_0.001,eqtls_non_0.001,eqtls_0.001),nrow=2,byrow=T)
    rows=c("non-eqtls","eqtls")
    cols=c("p<0.01","p>=0.01")
    dimnames(matrix) <- list("eQTL" = rows, "Selection" = cols)
    res0.001=oddsratio(matrix)
    
    eqtls_0.0001=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter((SNP %in% my_eqtls) & (SNP %in% snps0.0001)))$SNP))
    non_eqtls_0.0001=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter(!(SNP %in% my_eqtls) & (SNP %in% snps0.0001)))$SNP))
    eqtls_non_0.0001=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter((SNP %in% my_eqtls) & !(SNP %in% snps0.0001)))$SNP))
    non_eqtls_non_0.0001=length(unique((data_pruned_plus_HLA_eqtls_matchedAF %>% filter(!(SNP %in% my_eqtls) & !(SNP %in% snps0.0001)))$SNP))
    matrix=matrix(c(non_eqtls_non_0.0001,non_eqtls_0.0001,eqtls_non_0.0001,eqtls_0.0001),nrow=2,byrow=T)
    rows=c("non-eqtls","eqtls")
    cols=c("p<0.01","p>=0.01")
    dimnames(matrix) <- list("eQTL" = rows, "Selection" = cols)
    res0.0001=oddsratio(matrix)
    
    coeffs_groups=rbind(coeffs_groups,c(as.character(as.vector(annot)),
                                        res0.01$measure[2,1],res0.01$measure[2,2],res0.01$measure[2,3],
                                        res0.001$measure[2,1],res0.001$measure[2,2],res0.001$measure[2,3],
                                        res0.0001$measure[2,1],res0.0001$measure[2,2],res0.0001$measure[2,3]))
    Numbers=rbind(Numbers,c(eqtls_0.01,eqtls_0.001,eqtls_0.0001))
  }
  Altogether=cbind(coeffs_groups,Numbers)
  Altogether=data.frame(Altogether)
  colnames(Altogether)=c("Annotation","OR_1%","OR_0.025_1%","OR_0.975_1%","OR_0.1%","OR_0.025_0.1%","OR_0.975_0.1%",
                         "OR_0.01%","OR_0.025_0.01%","OR_0.975_0.01%","N_1%","N_0.1%","N_0.01%")
  
  for(i in 2:ncol(Altogether)){
    Altogether[,i]=as.numeric(as.vector(Altogether[,i]))
  }
  
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

### ENRICHMENT OF CIS-EQTLS (FROM QTLGEN, WHOLE BLOOD)
eqtls=fread("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/eQTL/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt")

#### Prune aDNA dataset to keep only independent variants (stronger on the HLA region)
AfterPrunning=as.character(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/GAIA/GENOMEWIDE2021/prunning100-10-0.6_maf0.01.prune.in",header=F,sep="\t")[,1]))
data_pruned=data %>% filter(SNP %in% AfterPrunning)
HLA_LDGroups=c("6:57","6:58","6:59","6:60","6:61","6:62","6:63","6:64","6:65")
stronger_pruned=as.character(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/GAIA/GENOMEWIDE2021/prunning1000-100-0.6_maf0.01.prune.in",header=F)[,1]))
ToInclude=(data %>% filter(GroupsLD %in% HLA_LDGroups & SNP %in% stronger_pruned))$SNP
data_pruned_tmp=data_pruned %>% filter(!(GroupsLD %in% HLA_LDGroups))
data_pruned_plus_HLA=rbind(data_pruned_tmp,data %>% filter(SNP %in% ToInclude))

### MATCH ON DAF
snps0.01_89loci=(data_pruned_plus_HLA %>% filter(pemp<0.01 & GroupsLD %in% groups))$SNP
groupsLD=(data_pruned_plus_HLA %>% filter(GroupsLD %in% groups))$SNP
BonferroniLevel=5e-08

eqtl_snps_89loci=(eqtls %>% filter(BonferroniP<BonferroniLevel & SNP %in% groupsLD))$SNP
Freq_distribution_eqtls_89loci=data_pruned_plus_HLA %>% filter(SNP %in% eqtl_snps_89loci) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
bins_eqtls_89loci=unique(Freq_distribution_eqtls_89loci$freq_bins)
TotalN_perbin_eqtls_89loci=(data_pruned_plus_HLA %>% filter(SNP %in% eqtl_snps_89loci) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(freq_bins) %>% summarise(NbrVars=length(pemp)) %>% filter(freq_bins %in% bins_eqtls_89loci))
TotalProp_perbin_eqtls_89loci=(data_pruned_plus_HLA %>% filter(SNP %in% eqtl_snps_89loci) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(freq_bins) %>% summarise(PropVars=length(pemp)/sum(TotalN_perbin_eqtls_89loci$NbrVars)) %>% filter(freq_bins %in% bins_eqtls_89loci))


eqtl_snps=(eqtls %>% filter(BonferroniP<BonferroniLevel))$SNP
Freq_distribution_eqtls=data_pruned_plus_HLA %>% filter(SNP %in% eqtl_snps) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
bins_eqtls=unique(Freq_distribution_eqtls$freq_bins)
TotalN_perbin_eqtls=(data_pruned_plus_HLA %>% filter(SNP %in% eqtl_snps) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(freq_bins) %>% summarise(NbrVars=length(pemp)) %>% filter(freq_bins %in% bins_eqtls))
TotalProp_perbin_eqtls=(data_pruned_plus_HLA %>% filter(SNP %in% eqtl_snps) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(freq_bins) %>% summarise(PropVars=length(pemp)/sum(TotalN_perbin_eqtls$NbrVars)) %>% filter(freq_bins %in% bins_eqtls))

Freq_distribution_allvars=data_pruned_plus_HLA %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
bins_allvars=Freq_distribution_allvars$freq_bins
TotalN_perbin_allvars=(data_pruned_plus_HLA %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(freq_bins) %>% summarise(NbrVars=length(pemp)) %>% filter(freq_bins %in% bins_eqtls))

## Find to maximum number of independent variants I could use to match the frequency distribution of the data I have
NbrNewVars=min(TotalN_perbin_allvars$NbrVars/TotalProp_perbin_eqtls$PropVars)
NbrNewVars_bin=floor(TotalProp_perbin_eqtls$PropVars*NbrNewVars)
set.seed(54)
count=1
vars_chosen=NULL
for(bin in as.character(as.vector(TotalN_perbin_allvars$freq_bins))){
  data_tmp=TotalN_perbin_allvars[which(TotalN_perbin_allvars$freq_bins==bin),]
  data_tmp_myvars=TotalN_perbin_eqtls[which(TotalN_perbin_eqtls$freq_bins==bin),]
  nbr_to_sample=NbrNewVars_bin[count]
  vars_chosen=c(vars_chosen,sample(Freq_distribution_allvars[which(Freq_distribution_allvars$freq_bins==bin),]$SNP,nbr_to_sample,replace=F))
  count=count+1
}
data_pruned_plus_HLA_eqtls_matchedAF=data_pruned_plus_HLA %>% filter(SNP %in% vars_chosen)

## function to compute enrichments
## first within candidate LD groups
snps0.01_89loci=(data_pruned_plus_HLA_eqtls_matchedAF %>% filter(pemp<0.01  & GroupsLD %in% groups))$SNP
snps0.001_89loci=(data_pruned_plus_HLA_eqtls_matchedAF %>% filter(pemp<0.001  & GroupsLD %in% groups))$SNP
snps0.0001_89loci=(data_pruned_plus_HLA_eqtls_matchedAF %>% filter(pemp<0.0001 & GroupsLD %in% groups))$SNP
within_groupsLD=enrichment_eqtl(snps0.01_89loci,snps0.001_89loci,snps0.0001_89loci)

## secondly, across all loci
snps0.01=(data_pruned_plus_HLA_eqtls_matchedAF %>% filter(pemp<0.01 ))$SNP
snps0.001=(data_pruned_plus_HLA_eqtls_matchedAF %>% filter(pemp<0.001 ))$SNP
snps0.0001=(data_pruned_plus_HLA_eqtls_matchedAF %>% filter(pemp<0.0001 ))$SNP
all_loci=enrichment_eqtl(snps0.01,snps0.001,snps0.0001)

## put together to plot
df=rbind(all_loci,within_groupsLD)
df$Data=c(rep("All loci",nrow(all_loci)),rep("89 loci",nrow(within_groupsLD)))
df$Data=factor(df$Data,levels=c("All loci","89 loci"))
df$labels=factor(df$labels,levels=rev(c("eQTL_8","eQTL_15","eQTL_50")))
df=data.frame(df)

## RESULT
df

