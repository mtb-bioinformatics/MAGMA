setwd("/Users/sandergoossens/Desktop/Project_Serial_samples")
library(readxl)
#read in VCF data coming from Lofreq ()
VCFData<-read_xls("/Users/sandergoossens/Desktop/Project_Serial_samples/B17W2testVCFInfo.xls")
# genTB data is imported after entering VCF from Lofreq from server to online into the tool and then taking the .var file.
# I made an excel from the GenTB.var output file doing the following in bash: wget url_of_outputfile
# then doing mv name_of_file.var name_of_files.xls
# name_of_file was in my case Sander.B017W2.dedup_reads.bam 
GenTBData<-read_xls("/Users/sandergoossens/Desktop/Project_Serial_samples/Sander.B017W2.dedup_reads.bam.xls")
#delete redundant info and merge files to get all info
AllData<-VCFData[,-c(3)]
AllData2<-GenTBData[,-c(4,5,9:15,19)]
AllData<-cbind(AllData,AllData2)
#Separate depth and AlleleFrequency
AllData[,"depth"]<-""
AllData[,"AF"]<-""
for (i in 1:length(AllData$CHROM)) {
  AllData$depth[i]<-sapply(strsplit(AllData$INFO[i], "DP="), "[", 2)
  AllData$depth[i]<-sapply(strsplit(AllData$depth[i], ";AF="), "[", 1)
  
  AllData$AF[i]<-sapply(strsplit(AllData$INFO[i], ";AF="), "[", 2)
  AllData$AF[i]<-sapply(strsplit(AllData$AF[i], ";SB="), "[", 1)
}
#load the list of coll with resistance mutations:
colllist<-read.delim("/Users/sandergoossens/Desktop/Project_Serial_samples/new_TBProfiler_panel_Coll2018(1).txt",header = TRUE, sep = "\t")
table(AllData$regionid1)
AllData[AllData$regionid1=="Rv0667",]
which(AllData$regionid1=="Rv0667")

#list of coll is incomplete
# so all mutations that GENTB tells me result in resistance should be callen in $resistance column
# due to low amount i did this manually atm, should code this better with lookup functions
AllData[,"Resistance"]<-"No"
# these mutations result in resistance, see json file GenTB: SNP_CN_2155168_C944G_S315T_katG en SNP_CN_761095_T1289C_L430P_rpoB;
# ik moet zien wat de andere mutaties zijn die in die lijst staan
# De andere rijen zijn de lineage SNPs, other en new snps die nog niet met resistance geassocieerd zijn, dus die kan ik erin laten
# Filter de important ones eruit
AllData$Resistance[AllData$varname=="SNP_CN_2155168_C944G_S315T_katG"]<-"INH"
AllData$Resistance[AllData$varname=="SNP_CN_761095_T1289C_L430P_rpoB"]<-"RIF"
length(table(AllData$desc1))

# Annotate pathways: step 1:  get the gene_IDs
# gene IDs are not in Entrez format
#gene_result.txt comes from GenTB--> file put in github as example
ListGeneIDs<-read.delim("/Users/sandergoossens/Desktop/gene_result.txt",header = TRUE, sep = "\t")
ListGeneIDs<-ListGeneIDs[,c(3,6:9,13:15)]
# write a convertor to add entrez_ID to AllData
AllData[,"EntrezGeneID"]<-""
ListGeneIDs$Aliases<-as.character(ListGeneIDs$Aliases)
# what to do with non-rv stuff, seems to be inter-Rv regions, these are not annotated, currently cut out?
# Make a code working only for AllData$regionid1 containing "rv".
table(AllData$regionid1[grepl("Rv",AllData$regionid1)])
#some genes are multiple times in the data --> DAVID will only use them as 1.
for (i in 1:length(AllData$CHROM)) {
  if (!is.na(match(AllData$regionid1[i], AllData$regionid1[grepl("Rv",AllData$regionid1)])) & length(which(grepl(AllData$regionid1[i],ListGeneIDs$Aliases)))>0) {
    AllData$EntrezGeneID[i]<-ListGeneIDs$GeneID[which(grepl(AllData$regionid1[i],ListGeneIDs$Aliases))]
  }
}
# make list to put in DAVID
AllData$EntrezGeneID<-as.numeric(AllData$EntrezGeneID)
which(is.na(AllData$EntrezGeneID))
listfortestingdavid<-AllData$EntrezGeneID[!is.na(AllData$EntrezGeneID)]
listfortestingdavid
clipr::write_clip((as.character(listfortestingdavid)))
# input for david is a list resulting from the code above
# gene_result.txt file with entrezIDs taken from NCBI gene as explained in: https://www.biostars.org/p/278673/
