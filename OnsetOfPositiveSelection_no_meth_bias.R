library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)

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

## create non-contiguous LD groups among the 89 loci by naming contiguous LD groups as only one of these
mydata = data
HLA_region=c("6:56","6:57","6:58","6:59","6:60","6:61","6:62","6:63","6:64")
mydata = mydata %>% mutate(GroupsLD = replace(GroupsLD,GroupsLD %in% c("19:39","19:40"), "19:39")) %>%
  mutate(GroupsLD = replace(GroupsLD,GroupsLD %in% c("2:184","2:185"), "2:184")) %>%
  mutate(GroupsLD = replace(GroupsLD,GroupsLD %in% HLA_region, "6:60")) %>%
  mutate(GroupsLD = replace(GroupsLD,GroupsLD %in% c("7:111","7:112"), "7:111")) %>%
  mutate(GroupsLD = replace(GroupsLD,GroupsLD %in% c("9:168","9:169"), "9:168")) %>%
  mutate(GroupsLD = replace(GroupsLD,GroupsLD %in% c("9:181","9:182"), "9:181"))


#### RANDOMLY SAMPLE SIMULATED VARIANTS MATCHED ON DAF AND SELECTION COEFFICIENT TO ASSESS POTENTIAL METHODOLOGICAL ARTIFACTS
## distributions of sel and DAF of the 89 top variants to match on with simulations
mydata2 = mydata %>% filter(pemp<0.01 & GroupsLD %in% groups)
mydata2_leadSNP = mydata2 %>%
  group_by(GroupsLD) %>%
  arrange(pbeta) %>%
  filter(row_number()==1)

ListLDGroups2=unique(mydata2_leadSNP$GroupsLD)
Onset=NULL
for(i in 1:length(ListLDGroups2)){
  Onset=c(Onset,(mydata2_leadSNP %>% filter(GroupsLD %in% ListLDGroups2[i]))$T)
}

Sel=NULL
for(i in 1:length(ListLDGroups2)){
  Sel=c(Sel,(mydata2_leadSNP %>% filter(GroupsLD %in% ListLDGroups2[i]))$s)
}

DAF=NULL
for(i in 1:length(ListLDGroups2)){
  DAF=c(DAF,(mydata2_leadSNP %>% filter(GroupsLD %in% ListLDGroups2[i]))$Epoch7.F)
}

## simulations matched on selection coefficient and DAF of the 89 top variants
SIMS=data.frame(fread("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ESTIMATION_WORKSPACE/ONSETOFSELECTION/tmp/RESmedianOnset_pvalmin.txt"))[1:100,]
colnames(SIMS)=paste("var",1:89,sep="")
SIMS$Groups=rep(1:100)
SIMS.m1 = reshape2::melt(SIMS, id.vars = c("Groups"),
                   measure.vars = colnames(SIMS)[1:(ncol(SIMS)-1)],variable.name="VarNbr",value.name="Medians")
Onset2=data.frame(Groups=rep(0,89),VarNbr=paste("var",1:89,sep=""),Medians=Onset)
SIMS.m1 = rbind(SIMS.m1,Onset2)
SIMS.m1=SIMS.m1 %>% mutate(Data=case_when(
  SIMS.m1$Groups %in% 1:100 ~ 'Simulations',
  TRUE ~ 'Real loci'))
SIMS.m1$Data <- factor(SIMS.m1$Data, levels = c("Simulations","Real loci"))

## simulations also matched on onset of selection of the 89 top variants
SIMS_matched=data.frame(fread("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ESTIMATION_WORKSPACE/ONSETOFSELECTION/tmp/RESmedianOnset_MatchedOnOnset_pvalmin.txt"))[1:100,]
colnames(SIMS_matched)=paste("var",1:89,sep="")
SIMS_matched$Groups=rep(101:200)
SIMS_matched.m1 = reshape2::melt(SIMS_matched, id.vars = c("Groups"),
                         measure.vars = colnames(SIMS_matched)[1:(ncol(SIMS_matched)-1)],variable.name="VarNbr",value.name="Medians")
SIMS_matched.m1$Data="Simulations_WithOnset"
SIMS_all=rbind(SIMS.m1,SIMS_matched.m1)
SIMS_all$Data <- factor(SIMS_all$Data, levels = c("Simulations","Simulations_WithOnset","Real loci"))

p=ggplot(SIMS_all,aes(x=Medians,group=Groups,linetype=Data,color=Data,size=Data))+
  geom_density(alpha=0.4)+
  scale_color_manual(values=c("grey50","grey80","black"))+
  scale_size_manual(values=c(0.5,0.5,2))+
  geom_density(data=SIMS_all %>% filter(Data=="Real loci"),aes(x=Medians,group=1))+
  # Custom the theme:
  theme_bw() +
  theme( 
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 15,angle=45,hjust=0.5,vjust=0.5),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    legend.key.width = unit(1.5,"cm"),
    legend.key.size = unit(1.5,"line")
  )


#### ONSET OF SELECTION COMPARISONS FOR ESTIMATIONS WITH AND WITHOUT MODERN DNA
mydata2_leadSNP_without_modern=read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/89_best/results_on_top_89_variants_no_modern_2.txt",header=F,sep="\t")
colnames(mydata2_leadSNP_without_modern)[2]="SNP"
colnames(mydata2_leadSNP_without_modern)[35]="T_nomodern"
mydata2_leadSNP_and_without_modern=merge(mydata2_leadSNP,mydata2_leadSNP_without_modern[,c("SNP","T_nomodern")],by="SNP")
plot(mydata2_leadSNP_and_without_modern$T,mydata2_leadSNP_and_without_modern$T_nomodern)
