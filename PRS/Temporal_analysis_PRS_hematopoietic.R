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
datas=loadRData(paste("~/Data/",FILE_NAME,sep=""))
groups=loadRData("~/Data/groupsLDEnrichment0.075pval0.01")
data=datas[[datatypes[4]]]
data=data %>% filter(!(abs(Epoch5.F-Epoch7.F)>0.1))


### LOAD ALSO NEGATIVE SELECTION
data_neg=datas[[datatypes[2]]]
data_neg=data_neg %>% filter(!(abs(Epoch7.F-Epoch5.F)>0.1))

#### TEMPORAL ANALYSIS OF HEMATOPOIETIC GWAS TRAITS

#### FIGURE 2D
### FIRST, CALCULATE PRS AND PVALUES FOR ALL 36 HEMATOPOIETIC TRAITS
### data was downloaded from here: https://www.phpc.cam.ac.uk/ceu/haematological-traits/

#### Prune aDNA dataset to keep only independent variants (stronger on the HLA region)
AfterPrunning=as.character(as.vector(read.table("~/Data/prunning100-10-0.6_maf0.01.prune.in",header=F,sep="\t")[,1]))
data_pruned=data %>% filter(SNP %in% AfterPrunning)
HLA_LDGroups=c("6:57","6:58","6:59","6:60","6:61","6:62","6:63","6:64","6:65")
stronger_pruned=as.character(as.vector(read.table("~/Data/prunning1000-100-0.6_maf0.01.prune.in",header=F)[,1]))
ToInclude=(data %>% filter(GroupsLD %in% HLA_LDGroups & SNP %in% stronger_pruned))$SNP
data_pruned_tmp=data_pruned %>% filter(!(GroupsLD %in% HLA_LDGroups))
data_pruned_plus_HLA=rbind(data_pruned_tmp,data %>% filter(SNP %in% ToInclude))

### LOAD ALL HEMATOPOIETIC-RELATED GWAS DATA (only genome-wide significant)
All_sig=fread("~/Data/ALL.gwas")
All_sig_vars=All_sig$VARIANT
All_sig=All_sig[-which(grepl("_",All_sig_vars)==F),]
All_sig_vars=All_sig_vars[-which(grepl("_",All_sig_vars)==F)]
chr=sapply(strsplit(All_sig_vars,split=":"),`[[`,1)
pos=sapply(strsplit(sapply(strsplit(All_sig_vars,split=":"),`[[`,2),split="_"),`[[`,1)
ref=sapply(strsplit(sapply(strsplit(All_sig_vars,split=":"),`[[`,2),split="_"),`[[`,2)
alt=sapply(strsplit(sapply(strsplit(All_sig_vars,split=":"),`[[`,2),split="_"),`[[`,3)
var1=paste(chr,pos,ref,alt,sep=":")
All=cbind(All_sig,var1)
colnames(All)[which(colnames(All)=="var1")]="CHR.POS.REF.ALT"
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 


##### NOW CALCULATE POLYGENIC SCORES AND ALL INFO FOR PLOTTING THE REST OF THE FIGURE

### FUNCTIONS
# to print r2 and pvalue on regression plots
lm_eqn <- function(df, y, x,z1,z2,z3,z4,z5,z6,wt){
  formula = as.formula(sprintf('%s ~ %s + %s+ %s+ %s+ %s+ %s+ %s', y, x,z1,z2,z3,z4,z5,z6))
  formula2 = as.formula(sprintf('%s ~ %s+ %s+ %s+ %s+ %s+ %s', y, z1,z2,z3,z4,z5,z6))
  m <- lm(formula, data=df, weights=wt);
  mfull=lm(formula, data=df, weights=Coverage)
  mnull=lm(formula2, data=df, weights=Coverage)
  # formating the values into a summary string to print out
  # ~ give some space, but equal size and comma need to be quoted
  eq <- substitute(~~italic(r)^2~"="~r2*","~~p~"="~italic(pvalue), 
                   list(target = y,
                        input = x,
                        a = format(as.vector(coef(m)[1]), digits = 2), 
                        b = format(as.vector(coef(m)[2]), digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3),
                        # getting the pvalue is painful
                        pvalue = format(anova(mnull, mfull, test = "Chisq")$`Pr(>Chi)`[2], digits=1)
                   )
  )
  as.character(as.expression(eq));                 
}



##### CALCULATE POLYGENIC SCORES FOR ALL TRAITS
fenos=names(table(All$baso))
names=names(table(All$baso))
names2=c("basophils","basophils_neutrophils_sum","basophils_p","basophils_p_granulocytes","eosinophils","eosinophils_basophils_sum","eosinophils_p",
         "eosinophils_p_granulocytes","granulocytes",
         "granulocytes_p_myeloid_white_cells","hematocrit","hemoglobin_concentration","high_light_scatter_reticulocyte","high_light_scatter_reticulocyte_p",
         "immature_fraction_reticulocytes","lymphocytes","lymphocytes_p","mean_corpuscular_hemoglobin","mean_corpuscular_hemoglobin_concentration",
         "mean_corpuscular_volume","monocytes","monocytes_p","mean_platelet_volume","myeloid_white_cell","neutrophils","neutrophils_eosinophils_sum",
         "neutrophils_p","neutrophils_p_granulocytes","plateletcrit","platelet_distribution_width","platelets","red_blood_cells","red_cell_distribution_width",
         "reticulocytes","reticulocytes_p","white_blood_cells")
FullFileName="~/Data/v44.3_1240K_public.anno.Extracted.Ancestries.txt.Extracted.Ancestries.txt"

### Initialize variables
pvaluesAll=NULL
pvaluesAfterBA=NULL
pvaluesBeforeBA=NULL
betaAll=NULL
betaBeforeBA=NULL
betaAfterBA=NULL
OUT_OUT=NULL
MATRIX=NULL
for(j in 1:length(fenos)){
  dataInds = read.csv(FullFileName,header=T,sep="\t")
  dataInds=data.frame(dataInds)
  NameAge="Date.mean.in.BP..OxCal.mu.for.a.direct.radiocarbon.date..and.average.of.range.for.a.contextual.date."
  NameID="Master.ID"
  NameCountry="Country"
  NameSource1="Anatolian"
  NameSource2="Yamnaya"
  NameSource3="Mesolithic_HG"
  # rename columns
  colnames(dataInds)[which(colnames(dataInds)==NameSource1)]="Anatolian"
  colnames(dataInds)[which(colnames(dataInds)==NameSource2)]="Yamnaya"
  colnames(dataInds)[which(colnames(dataInds)==NameAge)]="Age"
  colnames(dataInds)[which(colnames(dataInds)==NameID)]="ID"
  colnames(dataInds)[which(colnames(dataInds)==NameCountry)]="Country"
  
  ### HERE ELIMINATE ARMENIAN AND GEORGIAN INDIVIDUALS THAT REMAINED IF SO + RELATED INDIVIDUALS DEFINED BY READ OR METAFILE
  ListRelatedeness=read.table("~/Data/RELATEDENESS_BASED_ON_ANNOTATION_CAPTURE_1stDEGREE_INDS_TO_REMOVE.txt",header=F)
  rowsToEliminate=which(dataInds$Country=="Armenia" | dataInds$Country=="Georgia" | dataInds$ID %in% ListRelatedeness[,1])
  
  ## choose the hematopoeitic trait to study
  ra_hits=All %>% filter(baso==fenos[j])
  ra_hits=data.frame(ra_hits)
  ra_hits$`CHR.POS`=paste(ra_hits$CHR,ra_hits$BP,sep=":")
  
  ## work only in the 89 candidate loci
  data_reduced=data %>% filter(GroupsLD %in% groups)
  
  ## combine with positive selection data on the variants in aDNA (with a specific focus on the 89 loci)
  combined=merge(data_reduced,ra_hits,by="CHR.POS")
  LDGroups=unique(combined$GroupsLD)
  NbrHits=length(LDGroups)
  
  ## also combine with negative selection data to get both pvalues and s coefficients in case of need
  combined2=merge(combined,data_neg[,c("CHR.POS","pemp")],by="CHR.POS")
  
  ## here save the minimum psel value between positive and negative selection
  combined2$pemp=NA
  for(i in 1:nrow(combined2)){
    combined2$pemp[i]=min(combined2$pemp.x[i],combined2$pemp.y[i])
  }
  
  ## now format your data in a way that "OR.A1." is the OR of the allele that increases the value of the trait (so all OR>1)
  ## in A1 you have the effector allele (so the allele that increases the trait)
  combined2$OR.A1.=exp(as.numeric(as.vector(combined2$EFFECT)))
  tmp=combined2[which(combined2$OR.A1.<1),]$ALT.y
  combined2$A1=NA
  combined2$A2=NA
  
  ## if OR<1 that is because the effector allele (i.e., increasing the trait) should be the one noted as REF in GWAS
  ## if OR>1 that is because the effector allele (i.e., increasing the trait) should be the one noted as ALT in GWAS
  combined2[which(combined2$OR.A1.<1),]$A1=combined2[which(combined2$OR.A1.<1),]$REF.y
  combined2[which(combined2$OR.A1.<1),]$A2=combined2[which(combined2$OR.A1.<1),]$ALT.y
  combined2[which(combined2$OR.A1.>1),]$A1=combined2[which(combined2$OR.A1.>1),]$ALT.y
  combined2[which(combined2$OR.A1.>1),]$A2=combined2[which(combined2$OR.A1.>1),]$REF.y
  combined2$OR.original=combined2$OR.A1.
  combined2[which(combined2$OR.A1.<1),]$OR.A1.=1/combined2[which(combined2$OR.A1.<1),]$OR.A1.
  
  ## finally, frequencies in the aDNA data that usually are for derived alleles, I transform them here to the A1 allele
  combined2[which(combined2$A1==combined2$AA),c("Epoch2.F","Epoch3.F","Epoch4.F","Epoch5.F","Epoch6.F","Epoch7.F")]=1-combined2[which(combined2$A1==combined2$AA),c("Epoch2.F","Epoch3.F","Epoch4.F","Epoch5.F","Epoch6.F","Epoch7.F")]
  
  ## Here go one LD group by LD group and keep only the GWAS variant with the lowest GWAS p value
  hits=NULL
  for(i in 1:NbrHits){
    hits=rbind(hits,combined2[which(combined2$GroupsLD %in% LDGroups[i] & combined2$P==min(combined2[which(combined2$GroupsLD %in% LDGroups[i]),]$P)),][1,])
  }
  
  ## read data of all carriers of ped file for the allele that increases the trait
  ra_hits2=fread(paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/POLYGENIC/HEMA_TRAITS/hits_89loci_",names[j],"_individual_data.txt",sep=""),header=F)
  colnames_ra_hits2=as.character(as.vector(unlist(ra_hits2[1,])))
  colnames(ra_hits2)=as.character(as.vector(unlist(ra_hits2[1,])))
  ra_hits2=ra_hits2[-1,]
  ra_hits2=data.frame(ra_hits2)
  colnames_ra_hits2=sapply(strsplit(colnames_ra_hits2,split="_"),`[[`,1)
  
  ## keep columns to eliminate based on the criteria of only keeping europeans or non-first degree related individuals
  rows_to_eliminate_ra_hits2=which(colnames_ra_hits2 %in% as.character(as.vector(dataInds[rowsToEliminate,]$Version.ID)))

  ## in OUT save average frequencies for all variants to create a heatmap of the average frequency trajectory of the studied variants
  OUT=hits[1,]
  OUT$Epoch2.F=mean(hits$Epoch2.F);OUT$Epoch3.F=mean(hits$Epoch3.F);OUT$Epoch4.F=mean(hits$Epoch4.F);OUT$Epoch5.F=mean(hits$Epoch5.F);OUT$Epoch6.F=mean(hits$Epoch6.F);OUT$Epoch7.F=mean(hits$Epoch7.F)
  OUT$Gene=fenos[j]
  OUT$Epoch2.C=round(mean(hits$Epoch2.C));OUT$Epoch3.C=round(mean(hits$Epoch3.C));OUT$Epoch4.C=round(mean(hits$Epoch4.C));OUT$Epoch5.C=round(mean(hits$Epoch5.C));OUT$Epoch6.C=round(mean(hits$Epoch6.C));OUT$Epoch7.C=round(mean(hits$Epoch7.C))
  OUT_OUT=rbind(OUT_OUT,OUT)
  
  ## eliminate the individuals that are not europeans or that are first degree related
  colnames(ra_hits2)=colnames_ra_hits2
  dataInds=dataInds[-rowsToEliminate,]
  ra_hits2=ra_hits2[,-rows_to_eliminate_ra_hits2]
  colnames_ra_hits2=colnames_ra_hits2[-rows_to_eliminate_ra_hits2]
  
  ## also keep only individuals and no other info to create the Polygenic Scores
  colnames_ra_hits2=colnames_ra_hits2[-(1:6)]
  
  ## merge individual carrier data from the ped file with all the selection and GWAS data from before for each variant
  ra_hits2$`CHR.POS.COUNTED.ALT`=paste(ra_hits2$CHR,ra_hits2$POS,ra_hits2$COUNTED,ra_hits2$ALT,sep=":")
  merged=merge(ra_hits2,hits[,c("CHR.POS.COUNTED.ALT","A1","A2","OR.A1.","Epoch7.F")],by="CHR.POS.COUNTED.ALT")
  
  ## polarize the alleles so they all designate increasing trait (the counted allele in the file might not be the one increasing the trait)
  ## create a beta vector with effect sizes of each variant
  beta=NULL
  for(i in 1:nrow(merged)){
    merged[i,8:(ncol(merged)-4)]=ceiling(as.numeric(as.vector(merged[i,8:(ncol(merged)-4)]))/2)
    ## Either counted allele is the allele increasing the trait, then nothing to do but count carriers and average with beta
    if(merged[i,]$COUNTED==merged[i,]$A1){
      beta[i]=log(merged[i,]$OR.A1.)
      ## If that's not the case, then invert carriers so that the counted allele is the one increasing the trait
    }else{
      beta[i]=log(merged[i,]$OR.A1.)
      merged[i,8:(ncol(merged)-4)]=as.numeric(as.vector(merged[i,8:(ncol(merged)-4)]))+1
      merged[i,8:(ncol(merged)-4)][which(merged[i,8:(ncol(merged)-4)]==2)]=0
    }
  }
  
  # Calculate the PRS per individual as the sum of the effect sizes when carrier of "risk" allele divided by the total
  # sum of effect sizes in variants where the individual is covered
  RESPRS=apply(merged[,8:(ncol(merged)-4)],2,function(u) {sum(as.numeric(as.vector(u))*beta,na.rm=T)/sum(beta[which(!is.na(u))],na.rm=T)})
  
  ## keep also the coverage to then weight on it in the regression model
  COVPRS=apply(merged[,8:(ncol(merged)-4)],2,function(u) {length(which(!is.na(u)))})
  
  ## include the Version.ID with the computed scores to include individual info (Age, location and ancestry based on factor analysis)
  tmp=cbind(names(RESPRS),RESPRS,COVPRS)
  colnames(tmp)=c("Version.ID","Score","Coverage")
  tmp2=merge(tmp,dataInds[,c("Version.ID","Age","Lat.","Long.","PC1","PC2","PC3","PC4")],by="Version.ID")
  Toplot=tmp2
  Toplot=data.frame(Toplot)
  colnames(Toplot)=c("Subject.ID","Score","Coverage","Age","Lat","Long","F1","F2","F3","F4")
  Toplot$Score=as.numeric(as.vector(Toplot$Score))
  Toplot$Age=as.numeric(as.vector(Toplot$Age))
  Toplot$Coverage=as.numeric(as.vector(Toplot$Coverage))
  Toplot$Lat=as.numeric(as.vector(Toplot$Lat))
  Toplot$Long=as.numeric(as.vector(Toplot$Long))
  Toplot$F1=as.numeric(as.vector(Toplot$F1))
  Toplot$F2=as.numeric(as.vector(Toplot$F2))
  Toplot$F3=as.numeric(as.vector(Toplot$F3))
  Toplot$F4=as.numeric(as.vector(Toplot$F4))
  
  # Age put it as a negative value
  Toplot$Age=-Toplot$Age
  
  ## 3 categories
  ToplotAll=Toplot %>% filter(Age>(-10000))
  ToplotAfterBA=Toplot %>% filter(Age>(-4500))
  ToplotBeforeBA=Toplot %>% filter(Age>(-10000) & Age<(-4500))

  ## Save the p values for the model conducted on all 10,000 years or just after the BA
  pvaluesAll=rbind(pvaluesAll,c(names[j],lm_eqn(ToplotAll, 'Score','Age','Lat','Long','F1','F2',
                                                'F3','F4',ToplotAll$Coverage)))
  pvaluesBeforeBA=rbind(pvaluesBeforeBA,c(names[j],lm_eqn(ToplotBeforeBA, 'Score','Age','Lat','Long','F1','F2',
                                                          'F3','F4',ToplotBeforeBA$Coverage)))
  pvaluesAfterBA=rbind(pvaluesAfterBA,c(names[j],lm_eqn(ToplotAfterBA, 'Score','Age','Lat','Long','F1','F2',
                                                        'F3','F4',ToplotAfterBA$Coverage)))
  
  ## Idem beta values
  betaAll=c(betaAll,as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotAll, weights=ToplotAll$Coverage)$coefficients["Age"])))
  betaBeforeBA=c(betaBeforeBA,as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotBeforeBA, weights=ToplotBeforeBA$Coverage)$coefficients["Age"])))
  betaAfterBA=c(betaAfterBA,as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotAfterBA, weights=ToplotAfterBA$Coverage)$coefficients["Age"])))
  
  ## Average PRS per epoch
  PS1=mean((ToplotAll %>% filter(Age>(-10000) & Age<(-7500) & Coverage>0))$Score)
  PS2=mean((ToplotAll %>% filter(Age>(-7500) & Age<(-5000) & Coverage>0))$Score)
  PS3=mean((ToplotAll %>% filter(Age>(-5000) & Age<(-2500) & Coverage>0))$Score)
  PS4=mean((ToplotAll %>% filter(Age>(-2500) & Age<(0) & Coverage>0))$Score)
  
  ## in MATRIX save the heatmap on the four time transects for the PRS
  MATRIX=rbind(MATRIX,c(names2[j],PS1,PS2,PS3,PS4))
  cat(j,"\n")
}

MATRIX=data.frame(MATRIX)
colnames(MATRIX)=c("Trait","10-7.5 kya","7.5-5.0 kya","5.0-2.5 kya","2.5-0 kya")
MATRIX$Trait=names2

MATRIX_m1 = melt(MATRIX, id.vars = c("Trait"),
                               measure.vars = colnames(MATRIX)[2:ncol(MATRIX)])
colnames(MATRIX_m1)=c("Trait","Epoch","Score")
MATRIX_m1$Score=as.numeric(as.vector(MATRIX_m1$Score))
orden_pvals=order(as.numeric(as.vector(sapply(strsplit(pvaluesAll[,2],split="\""),`[[`,10))))
MATRIX_m1$Trait <- factor(MATRIX_m1$Trait, levels = rev(MATRIX_m1$Trait[orden_pvals]))
MATRIX_m1$Epoch <- factor(MATRIX_m1$Epoch, levels = c("10-7.5 kya","7.5-5.0 kya","5.0-2.5 kya","2.5-0 kya"))
## upper bound
MATRIX_m1$Score[which(MATRIX_m1$Score>0.66)]=0.66

### PLOTTING HEATMAP (ORDERED BY P VALUE ON THE LAST 10 MILLENNIA)
p_heat=ggplot(MATRIX_m1, aes(Epoch, Trait)) +
  geom_tile(aes(fill = Score), colour = "white") +
  scale_fill_gradient(limits=c(0.3,0.66),low = "white", high = "red")+
  ylab("")+
  theme_minimal() +
  theme( 
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size = 5),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30)
  )


##### hematological dot plot
hem=data.frame(Trait=names2,P=-log10(as.numeric(as.vector(sapply(strsplit(pvaluesAll[,2],split="\""),`[[`,10)))),
               P_BeforeBA=-log10(as.numeric(as.vector(sapply(strsplit(pvaluesBeforeBA[,2],split="\""),`[[`,10)))),
               P_AfterBA=-log10(as.numeric(as.vector(sapply(strsplit(pvaluesAfterBA[,2],split="\""),`[[`,10)))))
thres=-log10(0.05/(36))
hem$Trait <- factor(hem$Trait, levels = rev(hem$Trait[orden_pvals]))
colnames(hem)[which(colnames(hem)=="P")]="-log10(p)"

hem_m1 = melt(hem, id.vars = c("Trait"),
                 measure.vars = c("-log10(p)","P_BeforeBA","P_AfterBA"))
colnames(hem_m1)=c("Trait","Epoch","-log10(p)")
hem_m1$`-log10(p)`=as.numeric(as.vector(hem_m1$`-log10(p)`))
hem_m1$Epoch=as.character(as.vector(hem_m1$Epoch))

hem_m1=hem_m1 %>% mutate(Timing=case_when(
  grepl('log10',hem_m1$Epoch)~ 'All',
  grepl('BeforeBA',hem_m1$Epoch)~ 'Before_BA',
  grepl('AfterBA',hem_m1$Epoch)~ 'After_BA',
  TRUE ~ 'Other'),
  Significance=case_when(
    hem_m1$`-log10(p)`<thres~'Non-sig',
    TRUE ~ 'sig')
  )
hem_m1=hem_m1[rev(order(hem_m1$`-log10(p)`)),]

p_hem=ggplot(hem_m1, aes(x=Trait, y=`-log10(p)`, label=format(`-log10(p)`,digits=1))) + 
  geom_point(stat='identity', aes(col=Significance,fill=Significance,shape=Timing))  +
  scale_fill_manual(labels = c("Significant", "Non-significant"), 
                    values = c("sig"="darkblue", "Non-sig"="grey")) + 
  scale_color_manual(values = c("sig"="darkblue", "Non-sig"="grey")) +
  scale_shape_manual(values = c("All"=19, "Before_BA"=24,"After_BA"=25)) +
  theme_minimal()+
  ylab("")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 30),
        axis.text.x = element_text(size=20,angle=45,hjust=1,vjust=1),
        axis.text.y = element_blank(),
        plot.title = element_text(size=30),
        strip.text = element_text(size=30),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30))+
  coord_flip()+
  guides(fill = guide_legend(override.aes = list(size = 0)))
