### PRS OVER TIME FOR AUTOIMMUNITY

source("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/RSCRIPTS_PAPIER/Building_PRS_autoimmune_infection.R")

FullFileName="/Volumes/IGSR/GASPARD/DATING_SELECTION/GAIA/GENOMEWIDE2021/v44.3_1240K_public.anno.Extracted.Ancestries.txt.Extracted.Ancestries.txt"
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


#write.table(data_by_ldgroup_autoimmune[!duplicated(data_by_ldgroup_autoimmune$SNP),c("CHR.POS.REF.ALT","SNP")],paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/data_by_ldgroup_autoimmune.txt",sep=""),row.names=F,
#            col.names=T,quote=F,sep="\t")

#awk '{if(NR==FNR) {a[$1]=$1} else {if((($1":"$4":"$5":"$6 in a) || ($1":"$4":"$6":"$5 in a) || FNR==1)) {print $0}}}' /pasteur/zeus/projets/p02/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/data_by_ldgroup_autoimmune.txt /pasteur/zeus/projets/p02/IGSR/GASPARD/DATING_SELECTION/GAIA/GENOMEWIDE2021/v44.3_1240K_public.ann.vcf.gz.Extracted.traw > /pasteur/zeus/projets/p02/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/data_by_ldgroup_autoimmune_individual_data.txt

ra_hits2=fread(paste("/Volumes/IGSR/GASPARD/DATING_SELECTION/ZEUS/GWAS_CATALOG/data_by_ldgroup_autoimmune_individual_data.txt",sep=""),header=F)
colnames_ra_hits2=as.character(as.vector(unlist(ra_hits2[1,])))
colnames(ra_hits2)=as.character(as.vector(unlist(ra_hits2[1,])))
ra_hits2=ra_hits2[-1,]
ra_hits2=data.frame(ra_hits2)
colnames_ra_hits2=sapply(strsplit(colnames_ra_hits2,split="_"),`[[`,1)
rows_to_eliminate_ra_hits2=which(colnames_ra_hits2 %in% as.character(as.vector(dataInds[rowsToEliminate,]$Version.ID)))

colnames(ra_hits2)=colnames_ra_hits2

dataInds=dataInds[-rowsToEliminate,]
ra_hits2=ra_hits2[,-rows_to_eliminate_ra_hits2]
colnames_ra_hits2=colnames_ra_hits2[-rows_to_eliminate_ra_hits2]
colnames_ra_hits2=colnames_ra_hits2[-(1:6)]

ra_hits2$`CHR.POS.COUNTED.ALT`=paste(ra_hits2$CHR,ra_hits2$POS,ra_hits2$COUNTED,ra_hits2$ALT,sep=":")
merged=merge(ra_hits2,data_by_ldgroup_autoimmune[!duplicated(data_by_ldgroup_autoimmune$SNP),c("CHR.POS.COUNTED.ALT","A1","A2","OR.A1.","Epoch7.F")],by="CHR.POS.COUNTED.ALT")

## polarize the alleles so they all designate risk
## create a beta vector
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

RESPRS=apply(merged[,8:(ncol(merged)-4)],2,function(u) {sum(as.numeric(as.vector(u))*beta,na.rm=T)/sum(beta[which(!is.na(u))],na.rm=T)})
COVPRS=apply(merged[,8:(ncol(merged)-4)],2,function(u) {length(which(!is.na(u)))})

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


pvaluesAll=c("All autoimmune diseases",lm_eqn(ToplotAll, 'Score','Age','Lat','Long','F1','F2',
                                              'F3','F4',ToplotAll$Coverage))
pvaluesBeforeBA=c("All autoimmune diseases",lm_eqn(ToplotBeforeBA, 'Score','Age','Lat','Long','F1','F2',
                                                   'F3','F4',ToplotBeforeBA$Coverage))
pvaluesAfterBA=c("All autoimmune diseases",lm_eqn(ToplotAfterBA, 'Score','Age','Lat','Long','F1','F2',
                                                  'F3','F4',ToplotAfterBA$Coverage))

betaAll=as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotAll, weights=ToplotAll$Coverage)$coefficients["Age"]))
betaBeforeBA=as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotBeforeBA, weights=ToplotBeforeBA$Coverage)$coefficients["Age"]))
betaAfterBA=as.numeric(as.vector(lm(Score~Age+Lat+Long+F1+F2+F3+F4, data=ToplotAfterBA, weights=ToplotAfterBA$Coverage)$coefficients["Age"]))
