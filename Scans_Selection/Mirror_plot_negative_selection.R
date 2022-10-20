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

## LD groups found under our approach to be under positive selection
groups=loadRData("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ESTIMATION_WORKSPACE/tmp/groupsLDEnrichment0.075pval0.01")

## LD groups under positive selection to remove from the analysis of negative selection (including groups + classic loci under pos sel)
group_to_remove_neg=loadRData("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/ESTIMATION_WORKSPACE/tmp/group_to_remove_neg")

## Add GERP score to the variants
GERP=fread("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GERP/GERP_PositionsAnalysisNegSel.txt")
GERP=data.frame(GERP)
colnames(GERP)=c("CHR","POS","Neutral_score","gerp_RS")
GERP$`CHR.POS`=paste(GERP$CHR,GERP$POS,sep=":")
data_neg2=merge(data_neg,GERP[,c("CHR.POS","gerp_RS")],by="CHR.POS")
data_neg2=data_neg2[rev(order(data_neg2$s_l/data_neg2$thresholdSig0.0001)),]

## Show loci with at least 3 variants under neg sel
Significance=1e-02
NbrSigVariants=3
data_neg2_tmp = data_neg2 %>% filter(pemp<Significance)
SigLDGroups=sort(table(data_neg2_tmp$GroupsLD))
data_neg2_tmp=data_neg2_tmp %>% filter(GroupsLD %in% names(SigLDGroups[which(SigLDGroups>=NbrSigVariants)]))
ListLDGroups=unique(data_neg2_tmp$GroupsLD)
ListLDGroups=ListLDGroups[-which(ListLDGroups %in% group_to_remove_neg)]

## save candidate negatively-selected variants
conseq=(data_neg2 %>% filter(transcript_type=="protein_coding") )$Conseq
missense=names(table(conseq))[grepl("^missense",names(table(conseq)))]
snps_neg=(data_neg2 %>% filter( !(GroupsLD %in% group_to_remove_neg) & (Conseq %in% missense | Conseq=="structural_interaction_variant") & pemp<1e-02 & gerp_RS>4 ))$SNP
snps_neg_IIG=(data_neg2 %>% filter(!(GroupsLD %in% group_to_remove_neg) & (Conseq=="missense_variant" | Conseq=="structural_interaction_variant") & pemp<1e-02 & gerp_RS>4 & Genes %in% iig))$SNP


## MANHATTAN NEGATIVE
RelevantInfo=data_neg
RelevantInfo$RealLDGroups=RelevantInfo$GroupsLD
RelevantInfo2=RelevantInfo[,c("SNP","CHR","POS","pbeta","pemp","RealLDGroups","s","T")]
colnames(RelevantInfo2)[4]="P"
RelevantInfo2$P=as.numeric(as.vector(RelevantInfo2$P))

## !!!!!! ADDITION AFTER REVISIONS !!!!!! 
## correcting pvalues so to consider the minimum between the empirical p value and the pvalue from the theoretical approximation
RelevantInfo2 = RelevantInfo2 %>% mutate(pemp = ifelse(pemp==0, P, pemp))
RelevantInfo2$P=pmin(RelevantInfo2$P,RelevantInfo2$pemp)
## if pbeta is NA include pemp value
RelevantInfo2 = RelevantInfo2 %>% mutate(P = ifelse(is.na(P), pemp, P))
colnames(RelevantInfo2)[3]="BP"
RelevantInfo2$BP=as.numeric(as.vector(RelevantInfo2$BP))
##


don <- RelevantInfo2 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  #select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(RelevantInfo2, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse((RealLDGroups %in% ListLDGroups) | (SNP %in% snps_neg), "yes", "no")) %>%
  mutate( onset=cut(T, breaks=c(-Inf, 4500,Inf),labels=c("After BA","Before BA"))) %>%
  mutate( chromo=ifelse(CHR%%2==1,1,0)) %>%
  filter(-log10(pemp)>1)

onsets_by_LDGroups=don %>% filter(pemp<0.01) %>% group_by(RealLDGroups) %>% summarize(median_T=mean(T) )
new_don=merge(don,onsets_by_LDGroups,by="RealLDGroups")
axisdf = new_don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
new_don$median_T=as.numeric(as.vector(new_don$median_T))
new_don = new_don %>% mutate( mean_onset=cut(T, breaks=c(-Inf, 4500,Inf),labels=c("After BA","Before BA")))
new_don = new_don %>% mutate(candidate_IIG=ifelse(SNP %in% snps_neg_IIG,1,0))
new_don$candidate_IIG=as.factor(new_don$candidate_IIG)

## Use pastel colors for the top variant of each of the 89 loci
## palette pastel
library(scales)
red=brewer_pal(palette = "Set1")(8)[1]
yellow=brewer_pal(palette = "Pastel2")(8)[6]

p1=ggplot(subset(new_don %>% filter(-log10(P)>1), is_highlight=="no" ), aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), is_highlight=="no"  & (-log10(pemp))>(-log10(1e-04)) ), color="grey90", size=1.5) +
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), ((is_highlight=="no" & (-log10(pemp))<=(-log10(1e-04))) | (is_highlight=="yes")) & !(is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02))) & (chromo==1) ), col="#367BB4FF", alpha=0.8, size=2.5) +
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), ((is_highlight=="no" & (-log10(pemp))<=(-log10(1e-04))) | (is_highlight=="yes")) & !(is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02))) & (chromo==0) ), col="#0966B4FF", alpha=0.8, size=2.5) +
  
  
  
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), (is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02)) & (chromo==1)) ),  col="#367BB4FF",fill="#367BB4FF", alpha=0.8, size=2.5) +
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), (is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02)) & (chromo==0)) ),  col="#0966B4FF",fill="#0966B4FF", alpha=0.8, size=2.5) +
  
  
  geom_point(data=subset(new_don, (is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02)) & (SNP %in% snps_neg) & !(SNP %in% snps_neg_IIG)) ), aes(shape=mean_onset,col=mean_onset,fill=mean_onset), alpha=0.8, size=2.5) +
  geom_point(data=subset(new_don, (is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02)) & (SNP %in% snps_neg_IIG)) ), aes(shape=mean_onset,fill=mean_onset), col="black", alpha=0.8, size=2.5) +
  scale_color_manual(values=c(yellow,red)) +
  scale_shape_manual(values=c(25,24)) +
  scale_fill_manual(values=c(yellow,red)) +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ylim(c(0,20))+
  
  geom_hline(yintercept = 2,linetype = "longdash",col="black")+
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 15,angle=45,hjust=0.5,vjust=0.5),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

### s plot

p1_s=ggplot(subset(new_don %>% filter(-log10(P)>1), is_highlight=="no" ), aes(x=BPcum, y=-s)) +
  
  # Show all points
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), is_highlight=="no"  & (-log10(pemp))>(-log10(1e-04)) ), color="grey90", size=1.5) +
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), ((is_highlight=="no" & (-log10(pemp))<=(-log10(1e-04))) | (is_highlight=="yes")) & !(is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02))) & (chromo==1) ), col="#367BB4FF", alpha=0.8, size=2.5) +
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), ((is_highlight=="no" & (-log10(pemp))<=(-log10(1e-04))) | (is_highlight=="yes")) & !(is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02))) & (chromo==0) ), col="#0966B4FF", alpha=0.8, size=2.5) +
  
  
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), (is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02)) & (chromo==1)) ),  col="#367BB4FF",fill="#367BB4FF", alpha=0.8, size=2.5) +
  geom_point(data=subset(new_don%>% filter(-log10(P)>1), (is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02)) & (chromo==0)) ),  col="#0966B4FF",fill="#0966B4FF", alpha=0.8, size=2.5) +
  geom_point(data=subset(new_don, (is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02)) & (SNP %in% snps_neg) & !(SNP %in% snps_neg_IIG)) ), aes(shape=mean_onset,col=mean_onset,fill=mean_onset), alpha=0.8, size=2.5) +
  geom_point(data=subset(new_don, (is_highlight=="yes" & (-log10(pemp))>(-log10(1e-02)) & (SNP %in% snps_neg_IIG)) ), aes(shape=mean_onset,fill=mean_onset), col="black", alpha=0.8, size=2.5) +
  scale_color_manual(values=c(yellow,red)) +
  scale_shape_manual(values=c(25,24)) +
  scale_fill_manual(values=c(yellow,red)) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ylim(c(-0.1,0))+
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 15,angle=45,hjust=0.5,vjust=0.5),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )


#### END OF MIRROR PLOT
