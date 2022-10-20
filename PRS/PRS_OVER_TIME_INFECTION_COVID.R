### PRS OVER TIME FOR COVID
source("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/RSCRIPTS_PAPIER/Building_PRS_autoimmune_infection.R")

fenos=c("COVID_A2","COVID_B2")
names=fenos
names2=fenos
FullFileName="/Volumes/IGSR/GASPARD/DATING_SELECTION/GAIA/GENOMEWIDE2021/v44.3_1240K_public.anno.Extracted.Ancestries.txt.Extracted.Ancestries.txt"
pvaluesAll=NULL
pvaluesAfterBA=NULL
pvaluesBeforeBA=NULL
betaAll=NULL
betaBeforeBA=NULL
betaAfterBA=NULL
Nall=NULL
count=1
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
  
  
  ListRelatedeness=read.table("/Volumes/IGSR/GASPARD/DATING_SELECTION/GAIA/GENOMEWIDE2021/READ_RELATEDENESS/RELATEDENESS_BASED_ON_ANNOTATION_CAPTURE_1stDEGREE_INDS_TO_REMOVE.txt",header=F)
  rowsToEliminate=which(dataInds$Country=="Armenia" | dataInds$Country=="Georgia" | dataInds$ID %in% ListRelatedeness[,1])
  
  ## choose trait to study
  ra_hits=tian_and_covid %>% filter(phenotype==fenos[j])
  ra_hits=data.frame(ra_hits)
  
  
  ## now merge with aDNA
  combined=merge(data,ra_hits,by="SNP")
  combined$OR.A1.=combined$OR_risk
  combined$A1=combined$allele_risk
  combined$A2=apply(combined[,c("allele_risk","AA","DA")],1,function(i) {ifelse(i[1]==i[2],i[3],i[2])})
  ## once formatted, just extract the SNP with smalles GWAS p value for each LD group
  ## save info on negative selection if needed
  LDGroups=unique(combined$GroupsLD)
  NbrHits=length(LDGroups)
  combined2=merge(combined,data_neg[,c("SNP","pemp")],by="SNP")
  combined2$pemp=NA
  for(i in 1:nrow(combined2)){
    combined2$pemp[i]=min(combined2$pemp.x[i],combined2$pemp.y[i])
  }
  
  ## if the effector allele is the ancestral allele just update present allele frequency to the ancestral one 
  combined2[which(combined2$A1==combined2$AA),]$Epoch7.F=1-combined2[which(combined2$A1==combined2$AA),]$Epoch7.F
  
  hits=NULL
  for(i in 1:NbrHits){
    hits=rbind(hits,combined2[which(combined2$GroupsLD %in% LDGroups[i] & combined2$P==min(combined2[which(combined2$GroupsLD %in% LDGroups[i]),]$P)),][1,])
  }
  ## read data of all carriers of ped file for the allele that increases the risk
  ra_hits2=fread(paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/Tian/hits/hits_",names[j],"_individual_data.txt",sep=""),header=F)
  colnames_ra_hits2=as.character(as.vector(unlist(ra_hits2[1,])))
  colnames(ra_hits2)=as.character(as.vector(unlist(ra_hits2[1,])))
  ra_hits2=ra_hits2[-1,]
  ra_hits2=data.frame(ra_hits2)
  colnames_ra_hits2=sapply(strsplit(colnames_ra_hits2,split="_"),`[[`,1)
  rows_to_eliminate_ra_hits2=which(colnames_ra_hits2 %in% as.character(as.vector(dataInds[rowsToEliminate,]$Version.ID)))
  
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
  
  ## polarize the alleles so they all designate increasing risk (the counted allele in the file might not be the one increasing the risk)
  ## create a beta vector with effect sizes of each variant
  beta=NULL
  for(i in 1:nrow(merged)){
    ## First be sure that you got only pseudo-haploid data (0 or 1)
    merged[i,8:(ncol(merged)-4)]=ceiling(as.numeric(as.vector(merged[i,8:(ncol(merged)-4)]))/2)
    if(merged[i,]$COUNTED==merged[i,]$A1){
      beta[i]=log(merged[i,]$OR.A1.)
    }else{
      beta[i]=log(merged[i,]$OR.A1.)
      merged[i,8:(ncol(merged)-4)]=as.numeric(as.vector(merged[i,8:(ncol(merged)-4)]))+1
      merged[i,8:(ncol(merged)-4)][which(merged[i,8:(ncol(merged)-4)]==2)]=0
    }
  }
  # Calculate the PRS per individual as the sum of the effect sizes when carrier of risk allele divided by the total
  # sum of effect sizes in variants where the individual is covered
  
  #RESPRS=apply(merged[,8:(ncol(merged)-4)],2,function(u) {sum(as.numeric(as.vector(u))*beta,na.rm=T)/sum(beta[which(!is.na(u))])})
  RESPRS=apply(merged[,8:(ncol(merged)-4)],2,function(u) {sum(as.numeric(as.vector(u))*beta,na.rm=T)/sum(beta[which(!is.na(u))],na.rm=T)})
  
  ## keep also the coverage to then weight on it in the regression model
  COVPRS=apply(merged[,8:(ncol(merged)-4)],2,function(u) {length(which(!is.na(u)))})
  
  
  
  
  ## include the Version.ID with the computed scores to include individual info (Age, location and ancestry based on factor analysis)
  # by ancestry
  tmp=cbind(names(RESPRS),RESPRS,COVPRS)
  colnames(tmp)=c("Version.ID","Score","Coverage")
  tmp2=merge(tmp,dataInds[,c("Version.ID","Age","Lat.","Long.","PC1","PC2","PC3","PC4","Yamnaya","Anatolian","Mesolithic_HG","Country")],by="Version.ID")
  Toplot=tmp2
  Toplot=data.frame(Toplot)
  colnames(Toplot)=c("Subject.ID","Score","Coverage","Age","Lat","Long","F1","F2","F3","F4","Yamnaya","Anatolian","HG","Country")
  Toplot$Score=as.numeric(as.vector(Toplot$Score))
  Toplot$Age=as.numeric(as.vector(Toplot$Age))
  Toplot$Coverage=as.numeric(as.vector(Toplot$Coverage))
  Toplot$Lat=as.numeric(as.vector(Toplot$Lat))
  Toplot$Long=as.numeric(as.vector(Toplot$Long))
  Toplot$F1=as.numeric(as.vector(Toplot$F1))
  Toplot$F2=as.numeric(as.vector(Toplot$F2))
  Toplot$F3=as.numeric(as.vector(Toplot$F3))
  Toplot$F4=as.numeric(as.vector(Toplot$F4))
  Toplot$Yamnaya=as.numeric(as.vector(Toplot$Yamnaya))
  Toplot$Anatolian=as.numeric(as.vector(Toplot$Anatolian))
  Toplot$HG=as.numeric(as.vector(Toplot$HG))
  
  Toplot$Age=-Toplot$Age
  Toplottmp=Toplot
  
  Toplottmp=Toplottmp[-which(Toplottmp$Age<(-10000)),]
  Toplottmp$Source=NA
  Toplottmp[which(Toplottmp$HG>0.75),]$Source="HG"
  Toplottmp[which(Toplottmp$Anatolian>0.75 & Toplottmp$Age>(-10000)),]$Source="Anatolian"
  Toplottmp[which(Toplottmp$Yamnaya>0.75),]$Source="Yamnaya"
  
  Toplottmp[which(is.na(Toplottmp$Source)),]$Source="Admixed"
  Toplottmp=Toplottmp %>% filter(!is.na(Source) & !is.na(Score))
  Toplottmp$Source=factor(Toplottmp$Source, levels = c("Anatolian","HG","Yamnaya","Admixed"))
  
  
  ToplotAll =Toplottmp %>% filter(Age>(-10000))
  ToplotAfterBA =Toplottmp %>% filter(Age>(-4500))
  ToplotBeforeBA =Toplottmp %>% filter(Age<(-4500) & Age>(-10000))
  
  pvaluesAll=rbind(pvaluesAll,c(names2[j],lm_eqn(ToplotAll, 'Score','Age','Lat','Long','F1','F2',
                                                 'F3','F4',ToplotAll$Coverage)))
  pvaluesBeforeBA=rbind(pvaluesBeforeBA,c(names2[j],lm_eqn(ToplotBeforeBA, 'Score','Age','Lat','Long','F1','F2',
                                                           'F3','F4',ToplotBeforeBA$Coverage)))
  pvaluesAfterBA=rbind(pvaluesAfterBA,c(names2[j],lm_eqn(ToplotAfterBA, 'Score','Age','Lat','Long','F1','F2',
                                                         'F3','F4',ToplotAfterBA$Coverage)))
  
  betaAll=c(betaAll,as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotAll, weights=ToplotAll$Coverage)$coefficients["Age"])))
  betaBeforeBA=c(betaBeforeBA,as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotBeforeBA, weights=ToplotBeforeBA$Coverage)$coefficients["Age"])))
  betaAfterBA=c(betaAfterBA,as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotAfterBA, weights=ToplotAfterBA$Coverage)$coefficients["Age"])))
  
  Nall=c(Nall,nrow(hits))
}

### SUMMARIZE ALL INFO INTO A TABLE

INFECT=data.frame(Trait=names2,P=-log10(as.numeric(as.vector(sapply(strsplit(pvaluesAll[,2],split="\""),`[[`,10)))),
                  P_BeforeBA=-log10(as.numeric(as.vector(sapply(strsplit(pvaluesBeforeBA[,2],split="\""),`[[`,10)))),
                  P_AfterBA=-log10(as.numeric(as.vector(sapply(strsplit(pvaluesAfterBA[,2],split="\""),`[[`,10)))),
                  beta=betaAll,
                  beta_BeforeBA=betaBeforeBA,
                  beta_AfterBA=betaAfterBA,
                  N=Nall)
orden_pvals=order(as.numeric(as.vector(sapply(strsplit(pvaluesAll[,2],split="\""),`[[`,10))))
thres=-log10(0.05/(nrow(INFECT)))
INFECT$Trait <- factor(INFECT$Trait, levels = rev(INFECT$Trait[orden_pvals]))
colnames(INFECT)[which(colnames(INFECT)=="P")]="-log10(p)"
INFECT$Significance=ifelse(INFECT$`-log10(p)`<thres,"Non-sig","sig")
INFECT$Significance_BeforeBA=ifelse(INFECT$P_BeforeBA<thres,"Non-sig","sig")
INFECT$Significance_AfterBA=ifelse(INFECT$P_AfterBA<thres,"Non-sig","sig")
INFECT=INFECT[orden_pvals,]
INFECT2=rbind(INFECT,INFECT,INFECT)
INFECT2[(nrow(INFECT)+1):(nrow(INFECT)*2),]$`-log10(p)`=INFECT$P_BeforeBA
INFECT2[(nrow(INFECT)+1):(nrow(INFECT)*2),]$beta=INFECT$beta_BeforeBA
INFECT2[(nrow(INFECT)+1):(nrow(INFECT)*2),]$Significance=INFECT$Significance_BeforeBA
INFECT2[(nrow(INFECT)*2+1):(nrow(INFECT)*3),]$`-log10(p)`=INFECT$P_AfterBA
INFECT2[(nrow(INFECT)*2+1):(nrow(INFECT)*3),]$beta=INFECT$beta_AfterBA
INFECT2[(nrow(INFECT)*2+1):(nrow(INFECT)*3),]$Significance=INFECT$Significance_AfterBA
INFECT2$Timing=NA
INFECT2$Timing[1:nrow(INFECT)]="All"
INFECT2$Timing[(nrow(INFECT)+1):(nrow(INFECT)*2)]="Before_BA"
INFECT2$Timing[(nrow(INFECT)*2+1):(nrow(INFECT)*3)]="After_BA"
INFECT2$P_BeforeBA=NULL
INFECT2$P_AfterBA=NULL
INFECT2$beta_BeforeBA=NULL
INFECT2$beta_AfterBA=NULL
INFECT2$Significance_BeforeBA=NULL
INFECT2$Significance_AfterBA=NULL
INFECT2$p=10^(-INFECT2$`-log10(p)`)
