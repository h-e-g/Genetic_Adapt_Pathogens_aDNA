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
data_neg=datas[[datatypes[2]]]
colnames(data)[2]="CHR.POS"
colnames(data_neg)[2]="CHR.POS"
colnames(data)[3]="CHR.POS.COUNTED.ALT"
colnames(data_neg)[3]="CHR.POS.COUNTED.ALT"
data=data %>% filter(!(abs(Epoch5.F-Epoch7.F)>0.1))
data_neg=data_neg %>% filter(!(abs(Epoch7.F-Epoch5.F)>0.1))

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


### BUILD INFECTIOUS AND AUTOIMMUNE DATASETS
### RETRIEVE DATA FROM GWAS CATALOG AND TRY TO FIND SYSTEMATIC PLEIOTROPY BETWEEN INFECTIOUS AND AUTOIMMUNE DISORDERS
### GWAS catalog

## infectious traits
infectious_traits=read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/gwas-association-downloaded_2022-05-10-MONDO_0043544-withChildTraits_infections.tsv",header=T,sep="\t")
infectious_traits=data.frame(infectious_traits)

## autoimmune traits
autoimmune_traits=read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/gwas-association-downloaded_2022-05-10-MONDO_0005301-withChildTraits_autoimmune.tsv",header=T,sep="\t")
autoimmune_traits=data.frame(autoimmune_traits)

## keep only studies with European samples
infectious_traits=infectious_traits[grepl("Eur|eur",infectious_traits$INITIAL.SAMPLE.SIZE),]
autoimmune_traits=autoimmune_traits[grepl("Eur|eur",autoimmune_traits$INITIAL.SAMPLE.SIZE),]

## Only GWAS significant variants (Bonferroni)
infectious_traits$PVALUE_MLOG=as.numeric(as.vector(infectious_traits$PVALUE_MLOG))
infectious_traits=infectious_traits %>% filter(PVALUE_MLOG>(-log(5e-08,10)))
autoimmune_traits$PVALUE_MLOG=as.numeric(as.vector(autoimmune_traits$PVALUE_MLOG))
autoimmune_traits=autoimmune_traits %>% filter(PVALUE_MLOG>(-log(5e-08,10)))

## curate the list of traits
list_of_infections=as.character(as.vector(read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/infectious_traits_restricted.txt",header=F,sep="\t")[,1]))
list_of_autoimmune=as.character(as.vector(read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/autoimmune_traits_restricted.txt",header=F,sep="\t")[,1]))

### Restrict data to a list of curated traits
infectious_traits=infectious_traits %>% filter(DISEASE.TRAIT %in% list_of_infections)

### update OR variant from Patin et al., 2012
infectious_traits[which(infectious_traits$SNPS=="rs16851720"),]$OR.or.BETA=exp(0.39)

### If beta value provided instead of OR, transform it to an OR
infectious_traits$OR.or.BETA=as.numeric(as.vector(infectious_traits$OR.or.BETA))
infectious_traits[which(infectious_traits$OR.or.BETA<1 & grepl("unit",infectious_traits$X95..CI..TEXT.)),]$OR.or.BETA=exp(infectious_traits[which(infectious_traits$OR.or.BETA<1 & grepl("unit",infectious_traits$X95..CI..TEXT.)),]$OR.or.BETA)

### Then filter out variants without reported effect size
infectious_traits=infectious_traits %>% filter(!is.na(OR.or.BETA))

### And also filter out variants without reported risk allele
Variants_info=cbind(as.character(as.vector(infectious_traits$STRONGEST.SNP.RISK.ALLELE)),infectious_traits$PUBMEDID,as.character(as.vector(infectious_traits$DISEASE.TRAIT)))
Variants_info_rsID=sapply(strsplit(Variants_info[,1],split="-"),`[[`,1)
infectious_traits$risk_allele=sapply(strsplit(Variants_info[,1],split="-"),`[[`,2)
infectious_traits=infectious_traits %>% filter(!risk_allele=="?" & !grepl("[*]",infectious_traits$risk_allele))

### Restrict data to a list of curated traits
autoimmune_traits=autoimmune_traits %>% filter(DISEASE.TRAIT %in% list_of_autoimmune)

### If beta value provided instead of OR, transform it to an OR
autoimmune_traits$OR.or.BETA=as.numeric(as.vector(autoimmune_traits$OR.or.BETA))
autoimmune_traits[which(autoimmune_traits$OR.or.BETA<1 & grepl("unit",autoimmune_traits$X95..CI..TEXT.)),]$OR.or.BETA=exp(autoimmune_traits[which(autoimmune_traits$OR.or.BETA<1 & grepl("unit",autoimmune_traits$X95..CI..TEXT.)),]$OR.or.BETA)

### Then filter out variants without reported effect size
autoimmune_traits=autoimmune_traits %>% filter(!is.na(OR.or.BETA))

### And also filter out variants without reported risk allele
Variants_info_autoimmune=cbind(as.character(as.vector(autoimmune_traits$STRONGEST.SNP.RISK.ALLELE)),autoimmune_traits$PUBMEDID,as.character(as.vector(autoimmune_traits$DISEASE.TRAIT)))
Variants_info_autoimmune_rsID=sapply(strsplit(Variants_info_autoimmune[,1],split="-"),`[[`,1)
autoimmune_traits$risk_allele=sapply(strsplit(Variants_info_autoimmune[,1],split="-"),`[[`,2)
autoimmune_traits=autoimmune_traits %>% filter(!risk_allele=="?" & !grepl("[*]",autoimmune_traits$risk_allele))

### Combine infectious phenotypes with those from Tian et al. and the latest release of COVID, 
### and autoimmune phenotypes with those from CD, IBD and UC.
### Obtain lead SNPs for all infectious or autoimmune phenotypes together
### To that end, use 200 Kb windows to capture the lead SNP (in terms of GWAS p value of each window)

## COMBINE INFECTIOUS GWAS CATALOG PHENOTYPES WITH FULL GWAS DATA ON INFECTION FROM TIAN et al., AND THE COVID GWAS
## Tian
top8000=read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/Tian/Top8000Tian.csv",header=T,sep=";")
allele_noneff=sapply(strsplit(as.character(as.vector(top8000$alleles)),split="/"),`[[`,1)
allele_eff=sapply(strsplit(as.character(as.vector(top8000$alleles)),split="/"),`[[`,2)
top8000$allele_eff=allele_eff
top8000$allele_noneff=allele_noneff
top8000$OR=exp(top8000$effect)
top8000$allele_risk=ifelse(top8000$OR>1,allele_eff,allele_noneff)
top8000$allele_nonrisk=ifelse(top8000$OR>1,allele_noneff,allele_eff)
top8000$OR_risk=ifelse(top8000$OR>1,top8000$OR,1/top8000$OR)
top8000=data.frame(top8000)
top8000_bonferroni=top8000 %>% filter(pvalue<5e-08)
infectious_phenotypes=names(which(table(top8000_bonferroni$phenotype)>0))
colnames(top8000_bonferroni)[which(colnames(top8000_bonferroni)=="pvalue")]="P"

## COVID
covid_a2=read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/COVID19hg/COVID19_HGI_A2_ALL_20220403.10k_GRCh37_genome-wide_significant.tsv",header=T,sep="\t")
covid_a2=data.frame(covid_a2)
covid_a2$REF=as.character(as.vector(covid_a2$REF))
covid_a2$ALT=as.character(as.vector(covid_a2$ALT))
covid_a2$all_inv_var_meta_beta=as.numeric(as.vector(covid_a2$all_inv_var_meta_beta))
covid_a2$OR=exp(covid_a2$all_inv_var_meta_beta)
covid_a2$allele_risk=ifelse(covid_a2$OR>1,covid_a2$ALT,covid_a2$REF)
covid_a2$allele_nonrisk=ifelse(covid_a2$OR>1,covid_a2$REF,covid_a2$ALT)
covid_a2$OR_risk=ifelse(covid_a2$OR>1,covid_a2$OR,1/covid_a2$OR)
covid_a2$phenotype="COVID_A2"

covid_b1=read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/COVID19hg/COVID19_HGI_B1_ALL_20220403.10k_GRCh37_genome-wide_significant.tsv",header=T,sep="\t")
covid_b1=data.frame(covid_b1)
covid_b1$REF=as.character(as.vector(covid_b1$REF))
covid_b1$ALT=as.character(as.vector(covid_b1$ALT))
covid_b1$all_inv_var_meta_beta=as.numeric(as.vector(covid_b1$all_inv_var_meta_beta))
covid_b1$OR=exp(covid_b1$all_inv_var_meta_beta)
covid_b1$allele_risk=ifelse(covid_b1$OR>1,covid_b1$ALT,covid_b1$REF)
covid_b1$allele_nonrisk=ifelse(covid_b1$OR>1,covid_b1$REF,covid_b1$ALT)
covid_b1$OR_risk=ifelse(covid_b1$OR>1,covid_b1$OR,1/covid_b1$OR)
covid_b1$phenotype="COVID_B1"

covid_b2=read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/COVID19hg/COVID19_HGI_B2_ALL_20220403.10k_GRCh37_genome-wide_significant.tsv",header=T,sep="\t")
covid_b2=data.frame(covid_b2)
covid_b2$REF=as.character(as.vector(covid_b2$REF))
covid_b2$ALT=as.character(as.vector(covid_b2$ALT))
covid_b2$all_inv_var_meta_beta=as.numeric(as.vector(covid_b2$all_inv_var_meta_beta))
covid_b2$OR=exp(covid_b2$all_inv_var_meta_beta)
covid_b2$allele_risk=ifelse(covid_b2$OR>1,covid_b2$ALT,covid_b2$REF)
covid_b2$allele_nonrisk=ifelse(covid_b2$OR>1,covid_b2$REF,covid_b2$ALT)
covid_b2$OR_risk=ifelse(covid_b2$OR>1,covid_b2$OR,1/covid_b2$OR)
covid_b2$phenotype="COVID_B2"

covid_c2=read.csv("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/COVID19hg/COVID19_HGI_C2_ALL_20220403.10k_GRCh37_genome-wide_significant.tsv",header=T,sep="\t")
covid_c2=data.frame(covid_c2)
covid_c2$REF=as.character(as.vector(covid_c2$REF))
covid_c2$ALT=as.character(as.vector(covid_c2$ALT))
covid_c2$all_inv_var_meta_beta=as.numeric(as.vector(covid_c2$all_inv_var_meta_beta))
covid_c2$OR=exp(covid_c2$all_inv_var_meta_beta)
covid_c2$allele_risk=ifelse(covid_c2$OR>1,covid_c2$ALT,covid_c2$REF)
covid_c2$allele_nonrisk=ifelse(covid_c2$OR>1,covid_c2$REF,covid_c2$ALT)
covid_c2$OR_risk=ifelse(covid_c2$OR>1,covid_c2$OR,1/covid_c2$OR)
covid_c2$phenotype="COVID_C2"

covid=rbind(covid_a2,covid_b1,covid_b2,covid_c2)
colnames(covid)[which(colnames(covid)=="all_inv_var_meta_p")]="P"
covid=data.frame(covid)
infectious_phenotypes_tmp=c(unique(as.character(as.vector(covid$phenotype))),unique(as.character(as.vector(top8000_bonferroni$phenotype))))

## combine both groups
top8000_bonferroni_2=top8000_bonferroni[,c("assay.name","scaffold","position","allele_risk","OR_risk","phenotype","P")]
colnames(top8000_bonferroni_2)[1]="SNP"
top8000_bonferroni_2$scaffold=sapply(strsplit(as.character(as.vector(top8000_bonferroni_2$scaffold)),split="chr"),`[[`,2)
colnames(top8000_bonferroni_2)[2]="chr"
covid_2=covid[,c("rsid","X.CHR","POS","allele_risk","OR_risk","phenotype","P")]
colnames(covid_2)[2]="chr"
colnames(covid_2)[3]="position"
colnames(covid_2)[1]="SNP"

tian_and_covid=rbind(top8000_bonferroni_2,covid_2)
tian_and_covid=data.frame(tian_and_covid)
tian_and_covid$allele_risk=as.character(as.vector(tian_and_covid$allele_risk))
tian_and_covid$chr=as.numeric(as.vector(tian_and_covid$chr))
tian_and_covid$position=as.numeric(as.vector(tian_and_covid$position))
tian_and_covid$P=as.numeric(as.vector(tian_and_covid$P))
tian_and_covid=tian_and_covid[!grepl(":",tian_and_covid$SNP),]
tian_and_covid=tian_and_covid %>% filter(!is.na(SNP))

### HOMOGENEIZE WITH GWAS CATALOG BY TAKING LEAD SNPs PER PHENOTYPE (GWAS)
By_intervals_tian_and_covid=NULL
count=1
window_size=2e05
for(infect in unique(as.character(as.vector(tian_and_covid$phenotype)))){
  my_data_tmp=tian_and_covid %>% filter(phenotype==infect)
  if(any((my_data_tmp %>% group_by(chr) %>% summarise(largo=max(position)-min(position)))$largo==0)){
    chrs_with_single_variants=(my_data_tmp %>% group_by(chr) %>% summarise(largo=max(position)-min(position)) %>% filter(largo==0))$chr
    vars=as.character(as.vector((my_data_tmp %>% filter(chr %in% chrs_with_single_variants))$SNP))
    my_data_tmp=my_data_tmp %>% filter(!chr %in% chrs_with_single_variants)
  }else{
    vars=NULL
  }
  if(nrow(my_data_tmp)>0){
    vars=c(vars,unique(as.character(as.vector((my_data_tmp %>% group_by(chr) %>% mutate( intervals=cut(position, breaks=seq(min(position)-2e5,max(position)+2e5,min(window_size,max(position)-min(position))))) %>%
                                                 group_by(intervals) %>% filter(!is.na(intervals)) %>% arrange(P) %>% filter(row_number()==1))$SNP))))
  }
  
  By_intervals_tian_and_covid=rbind(By_intervals_tian_and_covid,my_data_tmp %>% filter(SNP %in% vars))
  cat(count,"\n")
  count=count+1
}
tian_and_covid_by_intervals=By_intervals_tian_and_covid

## get info on all infectious traits
infectious_traits_2=infectious_traits[,c("SNPS","CHR_ID","CHR_POS","STRONGEST.SNP.RISK.ALLELE","OR.or.BETA","DISEASE.TRAIT","PVALUE_MLOG")]
infectious_traits_2$STRONGEST.SNP.RISK.ALLELE=sapply(strsplit(as.character(as.vector(infectious_traits_2$STRONGEST.SNP.RISK.ALLELE)),split="-"),`[[`,2)
infectious_traits_2$PVALUE_MLOG=10^(-as.numeric(as.vector(infectious_traits_2$PVALUE_MLOG)))
colnames(infectious_traits_2)=colnames(tian_and_covid_by_intervals)
all_infectious_traits=rbind(tian_and_covid_by_intervals,infectious_traits_2)


## COMBINE AUTOIMMUNE GWAS CATALOG PHENOTYPES WITH FULL GWAS DATA ON AUTOIMMUNITY FROM GWAS ATLAS ON RA, UC, IBD AND CD
## get info on all autoimmune traits by combining IBD, UC and CD with the ones from GWAS catalog
fenos=c("RA_GWASmeta_European_v2_bonferroni.txt","cd_build37_40266_20161107_bonferroni.txt",
        "ibd_build37_59957_20161107_bonferroni.txt","uc_build37_45975_20161107_bonferroni.txt")

##### RA

ra_hits=data.frame(fread(paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/POLYGENIC/GWAS_ATLAS/",fenos[1],sep=""),header=T))
ra_hits$`CHR.POS`=paste(ra_hits$Chr,ra_hits$Position.hg19.,sep=":")
combined=merge(data,ra_hits,by="CHR.POS")
combined2=merge(combined,data_neg[,c("CHR.POS","pemp")],by="CHR.POS")
combined2$pemp=NA
for(i in 1:nrow(combined2)){
  combined2$pemp[i]=min(combined2$pemp.x[i],combined2$pemp.y[i])
}
tmp=combined2[which(combined2$OR.A1.<1),]$A1
combined2[which(combined2$OR.A1.<1),]$A1=combined2[which(combined2$OR.A1.<1),]$A2
combined2[which(combined2$OR.A1.<1),]$A2=tmp
combined2[which(combined2$OR.A1.<1),]$OR.A1.=1/combined2[which(combined2$OR.A1.<1),]$OR.A1.
combined2[which(combined2$A1==combined2$AA),]$Epoch7.F=1-combined2[which(combined2$A1==combined2$AA),]$Epoch7.F
combined2_ra=combined2

##### CD

cd_hits=data.frame(fread(paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/POLYGENIC/GWAS_ATLAS/",fenos[2],sep=""),header=T))
cd_hits$CHR=sapply(strsplit(as.character(as.vector(cd_hits$MarkerName)),split=":"),`[[`,1)
tmp=sapply(strsplit(as.character(as.vector(cd_hits$MarkerName)),split=":"),`[[`,2)
cd_hits$POS=sapply(strsplit(tmp,split="_"),`[[`,1)
cd_hits$REF=sapply(strsplit(tmp,split="_"),`[[`,2)
cd_hits$ALT=sapply(strsplit(tmp,split="_"),`[[`,3)
cd_hits$`CHR.POS.REF.ALT`=paste(cd_hits$CHR,cd_hits$POS,cd_hits$REF,cd_hits$ALT,sep=":")
combined=merge(data,cd_hits,by="CHR.POS.REF.ALT")
combined$Allele1[which(combined$Allele1=="a")]="A"
combined$Allele1[which(combined$Allele1=="c")]="C"
combined$Allele1[which(combined$Allele1=="g")]="G"
combined$Allele1[which(combined$Allele1=="t")]="T"
combined$Allele2[which(combined$Allele2=="a")]="A"
combined$Allele2[which(combined$Allele2=="c")]="C"
combined$Allele2[which(combined$Allele2=="g")]="G"
combined$Allele2[which(combined$Allele2=="t")]="T"
# A1 I define it as the risk allele
combined =combined %>% filter(Direction=="+++" | Direction=="---")
combined$A1=NA
combined$A2=NA
combined$OR.A1.=NA
combined[which(combined$Direction=="+++"),]$A1=combined[which(combined$Direction=="+++"),]$Allele2
combined[which(combined$Direction=="+++"),]$A2=combined[which(combined$Direction=="+++"),]$Allele1
combined[which(combined$Direction=="---"),]$A1=combined[which(combined$Direction=="---"),]$Allele1
combined[which(combined$Direction=="---"),]$A2=combined[which(combined$Direction=="---"),]$Allele2
combined[which(combined$Direction=="+++"),]$OR.A1.=exp(combined[which(combined$Direction=="+++"),]$Effect)
combined[which(combined$Direction=="---"),]$OR.A1.=1/exp(combined[which(combined$Direction=="---"),]$Effect)
combined2=merge(combined,data_neg[,c("CHR.POS.REF.ALT","pemp")],by="CHR.POS.REF.ALT")
combined2$pemp=NA
for(i in 1:nrow(combined2)){
  combined2$pemp[i]=min(combined2$pemp.x[i],combined2$pemp.y[i])
}
combined2[which(combined2$A1==combined2$AA),]$Epoch7.F=1-combined2[which(combined2$A1==combined2$AA),]$Epoch7.F
combined2_cd=combined2

##### IBD

ibd_hits=data.frame(fread(paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/POLYGENIC/GWAS_ATLAS/",fenos[3],sep=""),header=T))
ibd_hits$CHR=sapply(strsplit(as.character(as.vector(ibd_hits$MarkerName)),split=":"),`[[`,1)
tmp=sapply(strsplit(as.character(as.vector(ibd_hits$MarkerName)),split=":"),`[[`,2)
ibd_hits$POS=sapply(strsplit(tmp,split="_"),`[[`,1)
ibd_hits$REF=sapply(strsplit(tmp,split="_"),`[[`,2)
ibd_hits$ALT=sapply(strsplit(tmp,split="_"),`[[`,3)
ibd_hits$`CHR.POS.REF.ALT`=paste(ibd_hits$CHR,ibd_hits$POS,ibd_hits$REF,ibd_hits$ALT,sep=":")
combined=merge(data,ibd_hits,by="CHR.POS.REF.ALT")
combined$Allele1[which(combined$Allele1=="a")]="A"
combined$Allele1[which(combined$Allele1=="c")]="C"
combined$Allele1[which(combined$Allele1=="g")]="G"
combined$Allele1[which(combined$Allele1=="t")]="T"
combined$Allele2[which(combined$Allele2=="a")]="A"
combined$Allele2[which(combined$Allele2=="c")]="C"
combined$Allele2[which(combined$Allele2=="g")]="G"
combined$Allele2[which(combined$Allele2=="t")]="T"
# A1 I define it as the risk allele
combined =combined %>% filter(Direction=="+++" | Direction=="---")
combined$A1=NA
combined$A2=NA
combined$OR.A1.=NA
combined[which(combined$Direction=="+++"),]$A1=combined[which(combined$Direction=="+++"),]$Allele2
combined[which(combined$Direction=="+++"),]$A2=combined[which(combined$Direction=="+++"),]$Allele1
combined[which(combined$Direction=="---"),]$A1=combined[which(combined$Direction=="---"),]$Allele1
combined[which(combined$Direction=="---"),]$A2=combined[which(combined$Direction=="---"),]$Allele2
combined[which(combined$Direction=="+++"),]$OR.A1.=exp(combined[which(combined$Direction=="+++"),]$Effect)
combined[which(combined$Direction=="---"),]$OR.A1.=1/exp(combined[which(combined$Direction=="---"),]$Effect)
combined2=merge(combined,data_neg[,c("CHR.POS.REF.ALT","pemp")],by="CHR.POS.REF.ALT")
combined2$pemp=NA
for(i in 1:nrow(combined2)){
  combined2$pemp[i]=min(combined2$pemp.x[i],combined2$pemp.y[i])
}
combined2[which(combined2$A1==combined2$AA),]$Epoch7.F=1-combined2[which(combined2$A1==combined2$AA),]$Epoch7.F
combined2_ibd=combined2


##### UC

uc_hits=data.frame(fread(paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/POLYGENIC/GWAS_ATLAS/",fenos[4],sep=""),header=T))
uc_hits$CHR=sapply(strsplit(as.character(as.vector(uc_hits$MarkerName)),split=":"),`[[`,1)
tmp=sapply(strsplit(as.character(as.vector(uc_hits$MarkerName)),split=":"),`[[`,2)
uc_hits$POS=sapply(strsplit(tmp,split="_"),`[[`,1)
uc_hits$REF=sapply(strsplit(tmp,split="_"),`[[`,2)
uc_hits$ALT=sapply(strsplit(tmp,split="_"),`[[`,3)
uc_hits$`CHR.POS.REF.ALT`=paste(uc_hits$CHR,uc_hits$POS,uc_hits$REF,uc_hits$ALT,sep=":")
combined=merge(data,uc_hits,by="CHR.POS.REF.ALT")
combined$Allele1[which(combined$Allele1=="a")]="A"
combined$Allele1[which(combined$Allele1=="c")]="C"
combined$Allele1[which(combined$Allele1=="g")]="G"
combined$Allele1[which(combined$Allele1=="t")]="T"
combined$Allele2[which(combined$Allele2=="a")]="A"
combined$Allele2[which(combined$Allele2=="c")]="C"
combined$Allele2[which(combined$Allele2=="g")]="G"
combined$Allele2[which(combined$Allele2=="t")]="T"
# A1 I define it as the risk allele
combined =combined %>% filter(Direction=="+++" | Direction=="---")
combined$A1=NA
combined$A2=NA
combined$OR.A1.=NA
combined[which(combined$Direction=="+++"),]$A1=combined[which(combined$Direction=="+++"),]$Allele2
combined[which(combined$Direction=="+++"),]$A2=combined[which(combined$Direction=="+++"),]$Allele1
combined[which(combined$Direction=="---"),]$A1=combined[which(combined$Direction=="---"),]$Allele1
combined[which(combined$Direction=="---"),]$A2=combined[which(combined$Direction=="---"),]$Allele2
combined[which(combined$Direction=="+++"),]$OR.A1.=exp(combined[which(combined$Direction=="+++"),]$Effect)
combined[which(combined$Direction=="---"),]$OR.A1.=1/exp(combined[which(combined$Direction=="---"),]$Effect)
combined2=merge(combined,data_neg[,c("CHR.POS.REF.ALT","pemp")],by="CHR.POS.REF.ALT")
combined2$pemp=NA
for(i in 1:nrow(combined2)){
  combined2$pemp[i]=min(combined2$pemp.x[i],combined2$pemp.y[i])
}
combined2[which(combined2$A1==combined2$AA),]$Epoch7.F=1-combined2[which(combined2$A1==combined2$AA),]$Epoch7.F
combined2_uc=combined2

####

combined2_ra$phenotype="RA"
combined2_ibd$phenotype="IBD"
combined2_cd$phenotype="CD"
combined2_uc$phenotype="UC"

combined2_RA_2=combined2_ra[,c("CHR.POS.REF.ALT","SNP","A1","A2","OR.A1.")]
combined2_IBD_2=combined2_ibd[,c("CHR.POS.REF.ALT","SNP","A1","A2","OR.A1.")]
combined2_CD_2=combined2_cd[,c("CHR.POS.REF.ALT","SNP","A1","A2","OR.A1.")]
combined2_UC_2=combined2_uc[,c("CHR.POS.REF.ALT","SNP","A1","A2","OR.A1.")]

combined2_RA_2$pheno="RA"
combined2_IBD_2$pheno="IBD"
combined2_CD_2$pheno="CD"
combined2_UC_2$pheno="UC"
combined2_RA_CD_IBD_UC=rbind(combined2_RA_2,combined2_IBD_2,combined2_CD_2,combined2_UC_2)
SNPs_RA_CD_IBD_UC=c(combined2_RA_2$SNP,combined2_CD_2$SNP,combined2_UC_2$SNP,combined2_IBD_2$SNP)

tmp_RA=merge(combined2_RA_CD_IBD_UC[which(combined2_RA_CD_IBD_UC$pheno=="RA"),],combined2_ra[,c("SNP","P.val")],by="SNP")
colnames(tmp_RA)[which(colnames(tmp_RA)=="P.val")]="P.value"
tmp_IBD=merge(combined2_RA_CD_IBD_UC[which(combined2_RA_CD_IBD_UC$pheno=="IBD"),],combined2_ibd[,c("SNP","P.value")],by="SNP")
tmp_CD=merge(combined2_RA_CD_IBD_UC[which(combined2_RA_CD_IBD_UC$pheno=="CD"),],combined2_cd[,c("SNP","P.value")],by="SNP")
tmp_UC=merge(combined2_RA_CD_IBD_UC[which(combined2_RA_CD_IBD_UC$pheno=="UC"),],combined2_uc[,c("SNP","P.value")],by="SNP")

combined2_RA_CD_IBD_UC_with_p=rbind(tmp_RA,tmp_IBD,tmp_CD,tmp_UC)
colnames(combined2_RA_CD_IBD_UC_with_p)[length(colnames(combined2_RA_CD_IBD_UC_with_p))]="P"

combined2_RA_CD_IBD_UC_with_p$chr=sapply(strsplit(as.character(as.vector(combined2_RA_CD_IBD_UC_with_p$CHR.POS.REF.ALT)),split=":"),`[[`,1)
combined2_RA_CD_IBD_UC_with_p$position=sapply(strsplit(as.character(as.vector(combined2_RA_CD_IBD_UC_with_p$CHR.POS.REF.ALT)),split=":"),`[[`,2)
combined2_RA_CD_IBD_UC_with_p$chr=as.numeric(as.vector(combined2_RA_CD_IBD_UC_with_p$chr))
combined2_RA_CD_IBD_UC_with_p$position=as.numeric(as.vector(combined2_RA_CD_IBD_UC_with_p$position))
combined2_RA_CD_IBD_UC_with_p$P=as.numeric(as.vector(combined2_RA_CD_IBD_UC_with_p$P))


### HOMOGENEIZE WITH GWAS CATALOG BY TAKING LEAD SNPs PER PHENOTYPE (GWAS)

By_intervals_RA_CD_IBD_UC=NULL
count=1
window_size=2e05
for(inflam in unique(as.character(as.vector(combined2_RA_CD_IBD_UC_with_p$pheno)))){
  my_data_tmp=combined2_RA_CD_IBD_UC_with_p %>% filter(pheno==inflam)
  if(any((my_data_tmp %>% group_by(chr) %>% summarise(largo=max(position)-min(position)))$largo==0)){
    chrs_with_single_variants=(my_data_tmp %>% group_by(chr) %>% summarise(largo=max(position)-min(position)) %>% filter(largo==0))$chr
    vars=as.character(as.vector((my_data_tmp %>% filter(chr %in% chrs_with_single_variants))$SNP))
    my_data_tmp=my_data_tmp %>% filter(!chr %in% chrs_with_single_variants)
  }else{
    vars=NULL
  }
  if(nrow(my_data_tmp)>0){
    vars=c(vars,unique(as.character(as.vector((my_data_tmp %>% group_by(chr) %>% mutate( intervals=cut(position, breaks=seq(min(position)-2e5,max(position)+2e5,min(window_size,max(position)-min(position))))) %>%
                                                 group_by(intervals) %>% filter(!is.na(intervals)) %>% arrange(P) %>% filter(row_number()==1))$SNP))))
  }
  
  By_intervals_RA_CD_IBD_UC=rbind(By_intervals_RA_CD_IBD_UC,my_data_tmp %>% filter(SNP %in% vars))
  cat(count,"\n")
  count=count+1
}
RA_CD_IBD_UC_by_intervals=By_intervals_RA_CD_IBD_UC



autoimmune_traits_2=autoimmune_traits[,c("SNPS","CHR_ID","CHR_POS","STRONGEST.SNP.RISK.ALLELE","OR.or.BETA","DISEASE.TRAIT","PVALUE_MLOG")]
autoimmune_traits_2$STRONGEST.SNP.RISK.ALLELE=sapply(strsplit(as.character(as.vector(autoimmune_traits_2$STRONGEST.SNP.RISK.ALLELE)),split="-"),`[[`,2)
autoimmune_traits_2$PVALUE_MLOG=10^(-as.numeric(as.vector(autoimmune_traits_2$PVALUE_MLOG)))
colnames(autoimmune_traits_2)=colnames(tian_and_covid)
RA_CD_IBD_UC_by_intervals=RA_CD_IBD_UC_by_intervals[,c("SNP","chr","position","A1","OR.A1.","pheno","P")]
colnames(RA_CD_IBD_UC_by_intervals)=colnames(tian_and_covid)
all_autoimmune_traits=rbind(RA_CD_IBD_UC_by_intervals,autoimmune_traits_2)


## GET LEAD SNPS WHEN MERGING EITHER ALL AUTOIMMUNE OR ALL INFECTIOUS DISEASES TOGETHER (HERE USE EFFECT SIZE AND NOT P VALUE)
## THIS IS BECAUSE WE ARE COMBINING DIFFERENT GWAS WITH DIFFERENT SAMPLE SIZES

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# AUTOIMMUNE
By_intervals_all_autoimmune=NULL
count=1
window_size=2e05
all_autoimmune_traits$position=as.numeric(as.vector(all_autoimmune_traits$position))
my_data_tmp=all_autoimmune_traits %>% filter(!is.na(position) & !is.na(OR_risk) & chr %in% 1:22)

## Cases when for one chromosome I only have one variant
if(any((my_data_tmp %>% group_by(chr) %>% summarise(largo=max(position)-min(position)))$largo==0)){
  chrs_with_single_variants=(my_data_tmp %>% group_by(chr) %>% summarise(largo=max(position)-min(position)) %>% filter(largo==0))$chr
  vars=as.character(as.vector((my_data_tmp %>% filter(chr %in% chrs_with_single_variants))$SNP))
  my_data_tmp=my_data_tmp %>% filter(!chr %in% chrs_with_single_variants)
}else{
  vars=NULL
}
if(nrow(my_data_tmp)>0){
  #vars=c(vars,unique(as.character(as.vector((my_data_tmp %>% group_by(chr) %>% mutate( intervals=cut(position, breaks=seq(min(position)-2e5,max(position)+2e5,min(window_size,max(position)-min(position))))) %>%
  #                                             group_by(intervals) %>% filter(!is.na(intervals)) %>% arrange(P) %>% filter(row_number()==1))$SNP))))
  
  ## since I need to strongest effect size I need to put everything to the same scale
  ## So just for now I transform OR<1 to OR>1
  lines_OR_less_1=which(my_data_tmp$OR_risk<1)
  my_data_tmp[lines_OR_less_1,]$OR_risk=1/my_data_tmp[lines_OR_less_1,]$OR_risk
  vars=c(vars,unique(as.character(as.vector((my_data_tmp %>% group_by(chr) %>% mutate( intervals=cut(position, breaks=seq(min(position)-2e5,max(position)+2e5,min(window_size,max(position)-min(position))))) %>%
                                               group_by(intervals) %>% filter(!is.na(intervals)) %>% arrange(desc(OR_risk)) %>% filter(row_number()==1))$SNP))))
}

## Now I transform back the transformed ORs
my_data_tmp[lines_OR_less_1,]$OR_risk=1/my_data_tmp[lines_OR_less_1,]$OR_risk

## finally I keep the strongest (in terms of effect size) variant per window
By_intervals_all_autoimmune=my_data_tmp %>% filter(SNP %in% vars)
all_autoimmune_by_intervals=By_intervals_all_autoimmune

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Repeat for INFECTION
By_intervals_all_infectious=NULL
count=1
window_size=2e05
all_infectious_traits$position=as.numeric(as.vector(all_infectious_traits$position))
my_data_tmp=all_infectious_traits %>% filter(!is.na(position) & !is.na(OR_risk) & chr %in% 1:22)

## Cases when for one chromosome I only have one variant
if(any((my_data_tmp %>% group_by(chr) %>% summarise(largo=max(position)-min(position)))$largo==0)){
  chrs_with_single_variants=(my_data_tmp %>% group_by(chr) %>% summarise(largo=max(position)-min(position)) %>% filter(largo==0))$chr
  vars=as.character(as.vector((my_data_tmp %>% filter(chr %in% chrs_with_single_variants))$SNP))
  my_data_tmp=my_data_tmp %>% filter(!chr %in% chrs_with_single_variants)
}else{
  vars=NULL
}
if(nrow(my_data_tmp)>0){
  #vars=c(vars,unique(as.character(as.vector((my_data_tmp %>% group_by(chr) %>% mutate( intervals=cut(position, breaks=seq(min(position)-2e5,max(position)+2e5,min(window_size,max(position)-min(position))))) %>%
  #                                             group_by(intervals) %>% filter(!is.na(intervals)) %>% arrange(P) %>% filter(row_number()==1))$SNP))))
  
  ## since I need to strongest effect size I need to put everything to the same scale
  ## So just for now I transform OR<1 to OR>1
  lines_OR_less_1=which(my_data_tmp$OR_risk<1)
  my_data_tmp[lines_OR_less_1,]$OR_risk=1/my_data_tmp[lines_OR_less_1,]$OR_risk
  vars=c(vars,unique(as.character(as.vector((my_data_tmp %>% group_by(chr) %>% mutate( intervals=cut(position, breaks=seq(min(position)-2e5,max(position)+2e5,min(window_size,max(position)-min(position))))) %>%
                                               group_by(intervals) %>% filter(!is.na(intervals)) %>% arrange(desc(OR_risk)) %>% filter(row_number()==1))$SNP))))
  
}

## Now I transform back the transformed ORs
my_data_tmp[lines_OR_less_1,]$OR_risk=1/my_data_tmp[lines_OR_less_1,]$OR_risk

## finally I keep the strongest (in terms of effect size) variant per window
By_intervals_all_infectious=my_data_tmp %>% filter(SNP %in% vars)
all_infectious_by_intervals=By_intervals_all_infectious


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

##### UNCOMMENT TO TEST STATISTICAL EVIDENCES

# ## NOW USE ALL aDNA AS TAG SNPS, AND BUILD BLOCKS FOR TAG SNP USING WINDOWS OF 1MB AND R2 OF 0.6 AND MODERN DATA FROM EUR 1KG
# ## USING THESE BLOCKS WE WILL CHECK THE OVERLAP BETWEEN LEAD SNPS OF INFECTION AND AUTOIMMUNITY
# 
# tags_all=fread("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/snps_tags_0.6_1000_all.tags.list")
# tags_all=data.frame(tags_all)
# tags_all$TAGS2=paste(tags_all[,"SNP"],tags_all[,"TAGS"],sep="|")
# 
# ## get the blocks where you found at least one lead SNP for each sort of disorder
# haplos_infection=as.numeric(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/aDNA_haplotypes_including_infectious_GWAS_variants.txt",header=F)[,1]))
# haplos_autoimmune=as.numeric(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/aDNA_haplotypes_including_autoimmune_GWAS_variants.txt",header=F)[,1]))
# 
# tags_infection=tags_all[which(haplos_infection==1),]$SNP
# tags_autoimmune=tags_all[which(haplos_autoimmune==1),]$SNP
# tags_infection_autoimmune=intersect(tags_infection,tags_autoimmune)
# 
# ### Use plink to obtain pruned aDNA data by LD
# ### Calculate an OR based on pruned aDNA variants using plink (r2 of 0.6)
# ## r2 0.6 (only exclude variants if their r2>0.6)
# 
# INDEP_VARS=as.character(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/tags_aDNA_pruned_1000_0.6.prune.in",header=F)[,1]))
# tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6=intersect(unique(as.character(as.vector(tags_infection_autoimmune))),INDEP_VARS)
# tags_infection_GWAS_variants_aDNA_pruned_1000_0.6=intersect(unique(as.character(as.vector(tags_infection))),INDEP_VARS)
# tags_autoimmune_GWAS_variants_aDNA_pruned_1000_0.6=intersect(unique(as.character(as.vector(tags_autoimmune))),INDEP_VARS)
# tags_aDNA_pruned_1000_0.6=INDEP_VARS
# oddsratio(matrix(c(length(tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6),length(tags_infection_GWAS_variants_aDNA_pruned_1000_0.6)-length(tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6),
#                    length(tags_autoimmune_GWAS_variants_aDNA_pruned_1000_0.6)-length(tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6),length(tags_aDNA_pruned_1000_0.6)-length(tags_infection_GWAS_variants_aDNA_pruned_1000_0.6)-length(tags_autoimmune_GWAS_variants_aDNA_pruned_1000_0.6)+length(tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6)),nrow=2,byrow=T))
# 
# ## HERE WE SHOW THAT THESE PLEIOTROPIC VARIANTS ARE ENRICHED IN SELECTION
# ## WE USE RE-SAMPLING 
# bins=(data %>% filter(SNP %in% tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))))$freq_bins
# Freq_distribution=data %>% filter(SNP %in% INDEP_VARS) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
# Freq_distribution[which(Freq_distribution$Epoch7.F==0),"freq_bins"]="(0,0.05]"
# vars_chosen=NULL
# NbrRep=100
# set.seed(54)
# pvalues=NULL
# sel=NULL
# for(i in 1:NbrRep){
#   vars_chosen=NULL
#   for(j in 1:length(bins)){
#     vars_chosen=c(vars_chosen,sample(Freq_distribution[which(Freq_distribution$freq_bins==bins[j]),]$SNP,1))
#   }
#   merge_pos_neg_tmp=merge((data %>% filter(SNP %in% vars_chosen)),data_neg[,c("SNP","pemp","s")],by="SNP")
#   merge_pos_neg_tmp$pemp=apply(merge_pos_neg_tmp[,c("pemp.x","pemp.y")],1,min)
#   pvalues=c(pvalues,merge_pos_neg_tmp$pemp)
#   merge_pos_neg_tmp$s=apply(merge_pos_neg_tmp[,c("pemp.x","pemp.y","s.x","s.y")],1,function(i) {ifelse(i[1]<i[2],i[3],i[4])})
#   sel=c(sel,merge_pos_neg_tmp$s)
#   cat(i,"\n")
# }
# 
# merge_pos_neg=merge((data %>% filter(SNP %in% tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6)),data_neg[,c("SNP","pemp","s")],by="SNP")
# merge_pos_neg$pemp=apply(merge_pos_neg[,c("pemp.x","pemp.y")],1,min)
# merge_pos_neg$s=apply(merge_pos_neg[,c("pemp.x","pemp.y","s.x","s.y")],1,function(i) {ifelse(i[1]<i[2],i[3],i[4])})
# boxplot(pvalues,merge_pos_neg$pemp)
# wilcox.test(pvalues,merge_pos_neg$pemp,alternative = "greater")
# 
# 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 
# ### not allowing for the same snp to be tagged by more than one aDNA snp
# ## get the blocks where you found at least one lead SNP for each sort of disorder
# haplos_infection=as.numeric(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/aDNA_haplotypes_including_infectious_GWAS_variants2.txt",header=F)[,1]))
# haplos_autoimmune=as.numeric(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/aDNA_haplotypes_including_autoimmune_GWAS_variants2.txt",header=F)[,1]))
# 
# tags_infection=tags_all[which(haplos_infection==1),]$SNP
# tags_autoimmune=tags_all[which(haplos_autoimmune==1),]$SNP
# tags_infection_autoimmune=intersect(tags_infection,tags_autoimmune)
# 
# INDEP_VARS=as.character(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/tags_aDNA_pruned_1000_0.6.prune.in",header=F)[,1]))
# tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6=intersect(unique(as.character(as.vector(tags_infection_autoimmune))),INDEP_VARS)
# tags_infection_GWAS_variants_aDNA_pruned_1000_0.6=intersect(unique(as.character(as.vector(tags_infection))),INDEP_VARS)
# tags_autoimmune_GWAS_variants_aDNA_pruned_1000_0.6=intersect(unique(as.character(as.vector(tags_autoimmune))),INDEP_VARS)
# tags_aDNA_pruned_1000_0.6=INDEP_VARS
# oddsratio(matrix(c(length(tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6),length(tags_infection_GWAS_variants_aDNA_pruned_1000_0.6)-length(tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6),
#                    length(tags_autoimmune_GWAS_variants_aDNA_pruned_1000_0.6)-length(tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6),length(tags_aDNA_pruned_1000_0.6)-length(tags_infection_GWAS_variants_aDNA_pruned_1000_0.6)-length(tags_autoimmune_GWAS_variants_aDNA_pruned_1000_0.6)+length(tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6)),nrow=2,byrow=T))
# 
# ## HERE WE SHOW THAT THESE PLEIOTROPIC VARIANTS ARE ENRICHED IN SELECTION
# ## WE USE RE-SAMPLING 
# bins=(data %>% filter(SNP %in% tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))))$freq_bins
# Freq_distribution=data %>% filter(SNP %in% INDEP_VARS) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
# Freq_distribution[which(Freq_distribution$Epoch7.F==0),"freq_bins"]="(0,0.05]"
# vars_chosen=NULL
# NbrRep=1000
# set.seed(54)
# pvalues=NULL
# sel=NULL
# means_intersect=numeric(NbrRep)
# medians_intersect=numeric(NbrRep)
# for(i in 1:NbrRep){
#   vars_chosen=NULL
#   for(j in 1:length(bins)){
#     vars_chosen=c(vars_chosen,sample(Freq_distribution[which(Freq_distribution$freq_bins==bins[j]),]$SNP,1))
#   }
#   merge_pos_neg_tmp=merge((data %>% filter(SNP %in% vars_chosen)),data_neg[,c("SNP","pemp","s")],by="SNP")
#   merge_pos_neg_tmp$pemp=apply(merge_pos_neg_tmp[,c("pemp.x","pemp.y")],1,min)
#   pvalues=c(pvalues,merge_pos_neg_tmp$pemp)
#   merge_pos_neg_tmp$s=apply(merge_pos_neg_tmp[,c("pemp.x","pemp.y","s.x","s.y")],1,function(i) {ifelse(i[1]<i[2],i[3],i[4])})
#   sel=c(sel,merge_pos_neg_tmp$s)
#   means_intersect[i]=mean(merge_pos_neg_tmp$s)
#   medians_intersect[i]=median(merge_pos_neg_tmp$s)
#   cat(i,"\n")
# }
# 
# merge_pos_neg=merge((data %>% filter(SNP %in% tags_autoimmune_infection_GWAS_variants_aDNA_pruned_1000_0.6)),data_neg[,c("SNP","pemp","s")],by="SNP")
# merge_pos_neg$pemp=apply(merge_pos_neg[,c("pemp.x","pemp.y")],1,min)
# merge_pos_neg$s=apply(merge_pos_neg[,c("pemp.x","pemp.y","s.x","s.y")],1,function(i) {ifelse(i[1]<i[2],i[3],i[4])})
# sum(means_intersect>mean(merge_pos_neg$s))/NbrRep
# sum(medians_intersect>median(merge_pos_neg$s))/NbrRep
# ## empirical p value =  0 for medians and 0.019 for means after 1,000 reps
# wilcox.test(sel,merge_pos_neg$s,alternative = "less")
# 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 
# ### RESAMPLING OF VARIANTS MATCHED ON DAF AND # OF LD GROUPS
# ### now calculate significance of the overlap between infectious and autoimmune using resampling techniques matching on DAF and # of LD groups
# ## get groups of variants matching the structure of the tag variants for infection and autoimmunity
# # create bins of frequency and classify variants within these bins
# Freq_distribution=data %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
# Freq_distribution[which(Freq_distribution$Epoch7.F==0),"freq_bins"]="(0,0.05]"
# data_bins=data %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
# bins=unique(Freq_distribution$freq_bins)
# TotalN_perbin=(data %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(GroupsLD,freq_bins) %>% summarise(NbrVars=length(pemp)) %>% filter(freq_bins %in% bins))
# 
# ### For infection
# Freq_distribution_infection=data %>% filter(SNP %in% tags_infection) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
# Freq_distribution_infection[which(Freq_distribution_infection$Epoch7.F==0),"freq_bins"]="(0,0.05]"
# N_perbin_freq_infection=(data %>% filter(SNP %in% tags_infection) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(freq_bins) %>% summarise(NbrVars=length(pemp)))$NbrVars
# N_perbin_LD_infection=(data %>% filter(SNP %in% tags_infection) %>% group_by(GroupsLD) %>% summarise(NbrVars=length(pemp)))$NbrVars
# bins_infection=unique(Freq_distribution_infection$freq_bins)
# N_perbin_infection=(data %>% filter(SNP %in% tags_infection) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(GroupsLD,freq_bins) %>% summarise(NbrVars=length(pemp)) %>% filter(freq_bins %in% bins_infection))
# 
# ### For autoimmune
# Freq_distribution_autoimmune=data %>% filter(SNP %in% tags_autoimmune) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05)))
# Freq_distribution_autoimmune[which(Freq_distribution_autoimmune$Epoch7.F==0),"freq_bins"]="(0,0.05]"
# N_perbin_freq_autoimmune=(data %>% filter(SNP %in% tags_autoimmune) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(freq_bins) %>% summarise(NbrVars=length(pemp)))$NbrVars
# N_perbin_LD_autoimmune=(data %>% filter(SNP %in% tags_autoimmune) %>% group_by(GroupsLD) %>% summarise(NbrVars=length(pemp)))$NbrVars
# bins_autoimmune=Freq_distribution_autoimmune$freq_bins
# N_perbin_autoimmune=(data %>% filter(SNP %in% tags_autoimmune) %>% mutate(freq_bins=cut(Epoch7.F,breaks=seq(0,1,0.05))) %>% group_by(GroupsLD,freq_bins) %>% summarise(NbrVars=length(pemp)) %>% filter(freq_bins %in% bins_autoimmune))
# 
# ### RUN SCRIPT
# haplos_infection2=as.numeric(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/aDNA_haplotypes_including_infectious_GWAS_variants2.txt",header=F)[,1]))
# haplos_autoimmune2=as.numeric(as.vector(read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/DATA/aDNA_haplotypes_including_autoimmune_GWAS_variants2.txt",header=F)[,1]))
# 
# tags_infection2=tags_all[which(haplos_infection2==1),]$SNP
# tags_autoimmune2=tags_all[which(haplos_autoimmune2==1),]$SNP
# tags_infection_autoimmune2=intersect(tags_infection2,tags_autoimmune2)
# 
# infection_replicates=fread("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/infection_1000replicates2.txt",header=F,sep="\t")
# infection_replicates=data.frame(infection_replicates)
# autoimmune_replicates=fread("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/autoimmune_1000replicates2.txt",header=F,sep="\t")
# autoimmune_replicates=data.frame(autoimmune_replicates)
# 
# set.seed(54)
# NbrRep=1000
# length(intersect(tags_infection2,tags_autoimmune2))
# intersect_length=numeric(NbrRep)
# total_vars_infection=NULL
# total_vars_autoimmune=NULL
# pvalues_intersect=NULL
# sel_intersect=NULL
# means_intersect=numeric(NbrRep)
# medians_intersect=numeric(NbrRep)
# pvals_intersect=numeric(NbrRep)
# for(j in 1:NbrRep){
#   intersect_length[j]=length(intersect(as.character(as.vector(infection_replicates[j,])),as.character(as.vector(autoimmune_replicates[j,]))))
#   total_vars_infection=c(total_vars_infection,as.character(as.vector(infection_replicates[j,])))
#   total_vars_autoimmune=c(total_vars_autoimmune,as.character(as.vector(autoimmune_replicates[j,])))
#   
#   ## selection
#   merge_pos_neg_tmp=merge((data %>% filter(SNP %in% intersect(as.character(as.vector(infection_replicates[j,])),as.character(as.vector(autoimmune_replicates[j,]))))),data_neg[,c("SNP","pemp","s")],by="SNP")
#   merge_pos_neg_tmp$pemp=apply(merge_pos_neg_tmp[,c("pemp.x","pemp.y")],1,min)
#   pvalues_intersect=c(pvalues_intersect,merge_pos_neg_tmp$pemp)
#   merge_pos_neg_tmp$s=apply(merge_pos_neg_tmp[,c("pemp.x","pemp.y","s.x","s.y")],1,function(i) {ifelse(i[1]<i[2],i[3],i[4])})
#   sel_intersect=c(sel_intersect,merge_pos_neg_tmp$s)
#   means_intersect[j]=mean(merge_pos_neg_tmp$s)
#   medians_intersect[j]=median(merge_pos_neg_tmp$s)
#   pvals_intersect[j]=median(merge_pos_neg_tmp$pemp)
#   cat(j,"\n")
# }
# 
# # difference in terms of the ampount of overlap infection/autoimmune of variants between random replicates and real data
# merge_pos_neg=merge((data %>% filter(SNP %in% unique(as.character(as.vector(tags_infection_autoimmune2))))),data_neg[,c("SNP","pemp","s")],by="SNP")
# merge_pos_neg$pemp=apply(merge_pos_neg[,c("pemp.x","pemp.y")],1,min)
# merge_pos_neg$s=apply(merge_pos_neg[,c("pemp.x","pemp.y","s.x","s.y")],1,function(i) {ifelse(i[1]<i[2],i[3],i[4])})
# boxplot(pvalues_intersect,merge_pos_neg$pemp)
# boxplot(sel_intersect,merge_pos_neg$s)
# boxplot(medians_intersect,median(merge_pos_neg$s))
# boxplot(means_intersect,mean(merge_pos_neg$s))
# sum(means_intersect>mean(merge_pos_neg$s),na.rm=T)/length(which(!is.na(means_intersect)))
# sum(medians_intersect>median(merge_pos_neg$s),na.rm=T)/length(which(!is.na(medians_intersect)))
# 
# wilcox.test(sel_intersect,merge_pos_neg$s,alternative = "less")


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### BUILD PRS TRAJECTORY FOR ALL LEAD SNPS OF INFECTIOUS GWAS

ldgroups_infection=unique((data %>% filter(SNP %in% unique(as.character(as.vector(all_infectious_by_intervals$SNP)))))$GroupsLD)
all_data=merge(data,data_neg[,c("SNP","s","pemp")],by="SNP")
all_data$pemp=apply(all_data[,c("pemp.x","pemp.y")],1,min)
all_data$s=apply(all_data[,c("s.x","s.y","pemp.x","pemp.y")],1,function(i) {ifelse(i[3]<i[4],i[1],-i[2])})
data_infectious=merge(all_infectious_by_intervals,all_data[,c("SNP","pemp","GroupsLD","DA","s","CHR.POS.REF.ALT","CHR.POS.COUNTED.ALT","AA","Epoch7.F")],by="SNP")
data_by_ldgroup_infection=NULL
for(i in 1:length(ldgroups_infection)){
  #data_tmp=(data_infectious %>% filter(GroupsLD==ldgroups_infection[i]) %>% arrange(P))[1,]
  data_tmp=(data_infectious %>% filter(GroupsLD==ldgroups_infection[i]) %>% arrange(desc(OR_risk)))[1,]
  data_by_ldgroup_infection=rbind(data_by_ldgroup_infection,data_tmp)
}


## UPDATE RISK ALLELE OF rs16851720 ACCORDING TO PATIN ET AL., 2012
data_by_ldgroup_infection$A2=apply(data_by_ldgroup_infection[,c("allele_risk","DA","AA")],1,function(i) {ifelse(i[1]==i[2],i[3],i[2])})
data_by_ldgroup_infection$A1=data_by_ldgroup_infection$allele_risk
data_by_ldgroup_infection$OR.A1.=data_by_ldgroup_infection$OR_risk
data_by_ldgroup_infection[which(data_by_ldgroup_infection$SNP=="rs16851720"),c("OR.A1.")]=exp(0.39)
data_by_ldgroup_infection[which(data_by_ldgroup_infection$SNP=="rs16851720"),c("A1")]="A"
data_by_ldgroup_infection[which(data_by_ldgroup_infection$SNP=="rs16851720"),c("OR_risk")]=exp(0.39)
data_by_ldgroup_infection[which(data_by_ldgroup_infection$SNP=="rs16851720"),c("allele_risk")]="A"

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### BUILD PRS TRAJECTORY FOR ALL LEAD SNPS OF AUTO-IMMUNE GWAS MERGED TOGETHER
all_autoimmune_by_intervals2=all_autoimmune_by_intervals
ldgroups_autoimmune=unique((data %>% filter(SNP %in% unique(as.character(as.vector(all_autoimmune_by_intervals2$SNP)))))$GroupsLD)
all_data=merge(data,data_neg[,c("SNP","s","pemp")],by="SNP")
all_data$pemp=apply(all_data[,c("pemp.x","pemp.y")],1,min)
all_data$s=apply(all_data[,c("s.x","s.y","pemp.x","pemp.y")],1,function(i) {ifelse(i[3]<i[4],i[1],-i[2])})
data_autoimmune=merge(all_autoimmune_by_intervals2,all_data[,c("SNP","pemp","GroupsLD","DA","s","CHR.POS.REF.ALT","CHR.POS.COUNTED.ALT","AA","Epoch7.F")],by="SNP")
data_by_ldgroup_autoimmune=NULL
vars=NULL

lines_OR_less_1=which(data_autoimmune$OR_risk<1)
data_autoimmune[lines_OR_less_1,]$OR_risk=1/data_autoimmune[lines_OR_less_1,]$OR_risk
for(i in 1:length(ldgroups_autoimmune)){
  tmp=(data_autoimmune %>% filter(GroupsLD==ldgroups_autoimmune[i]) %>% arrange(desc(OR_risk)))[1,]
  vars=c(vars,tmp$SNP)
  #data_tmp=(tmp %>% filter(GroupsLD==ldgroups_autoimmune[i]) %>% arrange(desc(OR_risk)))[1,]
  #data_by_ldgroup_autoimmune=rbind(data_by_ldgroup_autoimmune,data_tmp)
}
data_autoimmune[lines_OR_less_1,]$OR_risk=1/data_autoimmune[lines_OR_less_1,]$OR_risk
data_by_ldgroup_autoimmune=data_autoimmune %>% filter(SNP %in% vars)

data_by_ldgroup_autoimmune$A2=apply(data_by_ldgroup_autoimmune[,c("allele_risk","DA","AA")],1,function(i) {ifelse(i[1]==i[2],i[3],i[2])})
data_by_ldgroup_autoimmune$A1=data_by_ldgroup_autoimmune$allele_risk
data_by_ldgroup_autoimmune$OR.A1.=data_by_ldgroup_autoimmune$OR_risk

## Change OR and allele_risk for those with OR_risk<1
data_by_ldgroup_autoimmune[which(data_by_ldgroup_autoimmune$OR_risk<1),]$A1=apply(data_by_ldgroup_autoimmune[which(data_by_ldgroup_autoimmune$OR_risk<1),c("allele_risk","DA","AA")],1,function(i) {ifelse(i[1]==i[2],i[3],i[2])})
data_by_ldgroup_autoimmune[which(data_by_ldgroup_autoimmune$OR_risk<1),]$A2=apply(data_by_ldgroup_autoimmune[which(data_by_ldgroup_autoimmune$OR_risk<1),c("allele_risk","DA","AA")],1,function(i) {ifelse(i[1]==i[2],i[2],i[3])})
data_by_ldgroup_autoimmune[which(data_by_ldgroup_autoimmune$OR_risk<1),]$OR.A1.=1/data_by_ldgroup_autoimmune[which(data_by_ldgroup_autoimmune$OR_risk<1),]$OR.A1.
