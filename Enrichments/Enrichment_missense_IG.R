library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggcorrplot)

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
data=data %>% filter(!(abs(Epoch5.F-Epoch7.F)>0.1))
groups=loadRData("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ESTIMATION_WORKSPACE/tmp/groupsLDEnrichment0.075pval0.01")

### THE DISTRIBUTION OF THE PROPORTION OF CANDIDATE VARIANTS FOR POSITIVE SELECTION IN LOCI INCLUDING 
### IIG WITH MISSENSE VARIANTS IS SIGNIFICANTLY HIGHER THAN THAT OF LOCI NOT FULFILLING THIS CONDITION
### THIS ANALYSIS FOCUSES ONLY ON THE 89 LOCI

### calculate enrichment of loci with positively-selected missense variants among loci with at least one iig (among all loci or the 89)

iig=as.character(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GO_ENRICHMENT/IIG.txt",header=F)[,1]))
miss_iig=unique((data %>% filter((Genes %in% iig) & Conseq %in% missense ))$GroupsLD)

Proportions=data %>% group_by(GroupsLD) %>% summarize(PropSig=sum(pemp<0.01)/length(pemp),Nbr=length(pemp),NbrSig=sum(pemp<0.01) )
TMP=Proportions %>% filter(Nbr>9 & NbrSig>2) %>% arrange(desc(PropSig))
wilcox.test((TMP %>% filter(GroupsLD %in% miss_iig & GroupsLD %in% groups))$PropSig,(TMP %>% filter(!(GroupsLD %in% miss_iig) & GroupsLD %in% groups))$PropSig,alternative = "greater")
