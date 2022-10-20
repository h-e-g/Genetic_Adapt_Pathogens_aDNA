### WORK AT THE VARIANT LEVEL TO GET VARIANTS SIGNIFICANTLY ASSOCIATED WITH THE AGE OF THE SAMPLES BY THEMSELVES (COVID)
source("~/PRS/Building_PRS_autoimmune_infection.R")

fenos=c("COVID_A2","COVID_B2")
names=c("COVID_A2","COVID_B2")
names2=c("COVID_A2","COVID_B2")
FullFileName="/Volumes/IGSR/GASPARD/DATING_SELECTION/GAIA/GENOMEWIDE2021/v44.3_1240K_public.anno.Extracted.Ancestries.txt.Extracted.Ancestries.txt"
BYVAR=list()
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
  
  BYVAR[[names2[j]]]=NULL
  for(var in 1:nrow(merged)){
    RESPRSbyvar=apply(merged[var,8:(ncol(merged)-4)],2,function(u) {sum(as.numeric(as.vector(u))*beta[var],na.rm=T)/sum(beta[var][which(!is.na(u))])})
    COVPRSbyvar=apply(merged[var,8:(ncol(merged)-4)],2,function(u) {length(which(!is.na(u)))})
    tmp=cbind(names(RESPRSbyvar),RESPRSbyvar,COVPRSbyvar)
    colnames(tmp)=c("Version.ID","Score","Coverage")
    tmp2=merge(tmp,dataInds[,c("Version.ID","Age","Anatolian","Yamnaya","Mesolithic_HG","Lat.","Long.","PC1","PC2","PC3","PC4")],by="Version.ID")
    Toplot=tmp2
    Toplot=data.frame(Toplot)
    colnames(Toplot)=c("Subject.ID","Score","Coverage","Age","Anatolian","Yamnaya","HG","Lat","Long","F1","F2","F3","F4")
    Toplot$Score=as.numeric(as.vector(Toplot$Score))
    Toplot$Age=as.numeric(as.vector(Toplot$Age))
    Toplot$Coverage=as.numeric(as.vector(Toplot$Coverage))
    Toplot$Anatolian=as.numeric(as.vector(Toplot$Anatolian))
    Toplot$Yamnaya=as.numeric(as.vector(Toplot$Yamnaya))
    Toplot$HG=as.numeric(as.vector(Toplot$HG))
    Toplot$Lat=as.numeric(as.vector(Toplot$Lat))
    Toplot$Long=as.numeric(as.vector(Toplot$Long))
    Toplot$F1=as.numeric(as.vector(Toplot$F1))
    Toplot$F2=as.numeric(as.vector(Toplot$F2))
    Toplot$F3=as.numeric(as.vector(Toplot$F3))
    Toplot$F4=as.numeric(as.vector(Toplot$F4))
    # eliminate older than 14,000
    Toplot=Toplot %>% filter(Coverage==1)
    Toplottmp=Toplot
    Toplottmp$Age=-Toplottmp$Age
    Toplottmp=Toplottmp[-which(Toplottmp$Age<(-14000)),]
    Toplottmp$Source=NA
    Toplottmp[which(Toplottmp$HG>0.75),]$Source="HG"
    Toplottmp[which(Toplottmp$Anatolian>0.75 & Toplottmp$Age>(-10000)),]$Source="Anatolian"
    Toplottmp[which(Toplottmp$Yamnaya>0.75),]$Source="Yamnaya"
    Toplottmp=Toplottmp %>% filter(!is.na(Source) & !is.na(Score))
    
    
    mfull=lm(Score~Age+F1+F2+F3+F4+Long+Lat, data=Toplottmp, weights=Coverage)
    mnull=lm(Score~F1+F2+F3+F4+Long+Lat, data=Toplottmp, weights=Coverage)
    pval=anova(mnull, mfull, test = "Chisq")[2,5]
    coeff=mfull$coeff["Age"]
    BYVAR[[names2[j]]]=rbind(BYVAR[[names2[j]]],c(merged[var,3],coeff,pval))
  }
  colnames(BYVAR[[names2[j]]])=c("SNP","coeff","pval")
  BYVAR[[names2[j]]]=data.frame(BYVAR[[names2[j]]])
  BYVAR[[names2[j]]]$pval=as.numeric(as.vector(BYVAR[[names2[j]]]$pval))
  cat(j,"\n")
}

TotalNbrVars=0
for(j in 1:length(names2)){
  TotalNbrVars=TotalNbrVars+nrow(BYVAR[[names2[j]]])
}

## Get all significant variants after multiple testing correction
SigVars=NULL
for(j in 1:length(names2)){
  SigVars=rbind(SigVars,BYVAR[[names2[j]]] %>% filter(pval<(0.05/(TotalNbrVars))))
}

# SNP                 coeff         pval
# 1  rs6489864 -3.45488679947426e-05 0.0004161149
# 2 rs10774679 -3.61250600093044e-05 0.0002052213
