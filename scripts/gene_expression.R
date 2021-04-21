


# author:
# date:
# about:


#setwd('/extData/NGS/hurlab/temi/projects/camda/reproducible/scripts/')


curr <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(curr))

# load packages
source('../scripts/packages.R')

data_dir <- '../data/'

##Load in Data
#Drug Properties
mold2 <-read.csv("../data/p1-mold2-camda2020.csv")
faers <-read.csv("../data/p2-faers-camda2020.csv")
tox21 <-read.csv("../data/p9-tox21-camda2020.csv")

#L1000 Z score expression Data from Cell testing
phh <-read.csv("../data/p3-phh-camda2020.csv")
hepg2 <-read.csv("../data/p4-hepg2-camda2020.csv")
ha1e <-read.csv("../data/p5-haie-camda2020.csv")
a375 <-read.csv("../data/p6-a375-camda2020.csv")
mcf7 <-read.csv("../data/p7-mcf7-camda2020.csv")
pc3 <-read.csv("../data/p8-pc3-camda2020.csv")

#Negative Ctl
random <- read.csv("../data/p10-random-camda2020.csv")

#Training and Validation
targets <- read.csv("../data/targets-camda2020.csv")
Training_CMID <- targets[which(targets$Training_Validation=="Training Set"),]

##############
phh3 <- phh
row.names(phh3) <- phh$CAM_ID
phh3 <- phh3[,-1]

hepg2_3 <- hepg2
row.names(hepg2_3) <- hepg2$CAM_ID
hepg2_3 <- hepg2_3[,-1]

ha1e3 <- ha1e
row.names(ha1e3) <- ha1e$CAM_ID
ha1e3 <- ha1e3[,-1]

a375_3 <- a375
row.names(a375_3) <- a375$CAM_ID
a375_3 <- a375_3[,-1]

mcf7_3 <- mcf7
row.names(mcf7_3) <- mcf7$CAM_ID
mcf7_3 <- mcf7_3[,-1]

pc3_3 <- pc3
row.names(pc3_3) <- pc3$CAM_ID
pc3_3 <- pc3_3[,-1]
phh_rankx <- as.data.frame(apply(phh3, 1, rank, ties.method="average"))

#add _phh to each drug
phh_rankx <- phh_rankx %>% dplyr::rename_all(paste0,"_phh")
phh_rankx2 <- as.data.frame(t(phh_rankx))

#Add sample info
phh_rankx2$instance_id <- rownames(phh_rankx2)
phh_rankx2$cell <- "PHH"
phh_rankx2$normalized_name <- rownames(phh3)

phh_sample_info <- as.data.frame(phh_rankx2$instance_id)
colnames(phh_sample_info) <- "instance_id"
phh_sample_info$normalized_name <- phh_rankx2$normalized_name
phh_sample_info$cell <- phh_rankx2$cell


##
hepg2_3 <- hepg2
row.names(hepg2_3) <- hepg2$CAM_ID
hepg2_3 <- hepg2_3[,-1]
hepg2_rankx <- as.data.frame(apply(hepg2_3, 1, rank, ties.method="average"))

hepg2_rankx <- hepg2_rankx %>% dplyr::rename_all(paste0, "_hepg2")
hepg2_rankx2 <- as.data.frame(t(hepg2_rankx))

hepg2_rankx2$instance_id <- rownames(hepg2_rankx2)
hepg2_rankx2$cell <- "HEPG2"
hepg2_rankx2$normalized_name <- rownames(hepg2_3)

hepg2_sample_info <- as.data.frame(hepg2_rankx2$instance_id)
colnames(hepg2_sample_info) <- "instance_id"
hepg2_sample_info$normalized_name <- hepg2_rankx2$normalized_name
hepg2_sample_info$cell <- hepg2_rankx2$cell

##
ha1e3 <- ha1e
row.names(ha1e3) <- ha1e$CAM_ID
ha1e3 <- ha1e3[,-1]
ha1e_rankx <- as.data.frame(apply(ha1e3, 1, rank, ties.method="average"))

ha1e_rankx <- ha1e_rankx %>% dplyr::rename_all(paste0, "_ha1e")
ha1e_rankx2 <- as.data.frame(t(ha1e_rankx))

ha1e_rankx2$instance_id <- rownames(ha1e_rankx2)
ha1e_rankx2$cell <- "HA1E"
ha1e_rankx2$normalized_name <- rownames(ha1e3)

ha1e_sample_info <- as.data.frame(ha1e_rankx2$instance_id)
colnames(ha1e_sample_info) <- "instance_id"
ha1e_sample_info$normalized_name <- ha1e_rankx2$normalized_name
ha1e_sample_info$cell <- ha1e_rankx2$cell


##
a375_3 <- a375
row.names(a375_3) <- a375$CAM_ID
a375_3 <- a375_3[,-1]
a375_rankx <- as.data.frame(apply(a375_3, 1, rank, ties.method="average"))

a375_rankx <- a375_rankx %>% dplyr::rename_all(paste0, "_a375")
a375_rankx2 <- as.data.frame(t(a375_rankx))

a375_rankx2$instance_id <- rownames(a375_rankx2)
a375_rankx2$cell <- "A375"
a375_rankx2$normalized_name <- rownames(a375_3)

a375_sample_info <- as.data.frame(a375_rankx2$instance_id)
colnames(a375_sample_info) <- "instance_id"
a375_sample_info$normalized_name <- a375_rankx2$normalized_name
a375_sample_info$cell <- a375_rankx2$cell


##
mcf7_3 <- mcf7
row.names(mcf7_3) <- mcf7$CAM_ID
mcf7_3 <- mcf7_3[,-1]
mcf7_rankx <- as.data.frame(apply(mcf7_3, 1, rank, ties.method="average"))

mcf7_rankx <- mcf7_rankx %>% dplyr::rename_all(paste0, "_mcf7")

mcf7_rankx2 <- as.data.frame(t(mcf7_rankx))

mcf7_rankx2$instance_id <- rownames(mcf7_rankx2)
mcf7_rankx2$cell <- "MCF7"
mcf7_rankx2$normalized_name <- rownames(mcf7_3)

mcf7_sample_info <- as.data.frame(mcf7_rankx2$instance_id)
colnames(mcf7_sample_info) <- "instance_id"
mcf7_sample_info$normalized_name<- mcf7_rankx2$normalized_name
mcf7_sample_info$cell <- mcf7_rankx2$cell

##
pc3_3 <- pc3
row.names(pc3_3) <- pc3$CAM_ID
pc3_3 <- pc3_3[,-1]
pc3_rankx <- as.data.frame(apply(pc3_3, 1, rank, ties.method="average"))

pc3_rankx <- pc3_rankx %>% dplyr::rename_all(paste0, "_pc3")
pc3_rankx2 <- as.data.frame(t(pc3_rankx))

pc3_rankx2$instance_id <- rownames(pc3_rankx2)
pc3_rankx2$cell <- "PC3"
pc3_rankx2$normalized_name <- rownames(pc3_3)

pc3_sample_info <- as.data.frame(pc3_rankx2$instance_id)
colnames(pc3_sample_info) <- "instance_id"
pc3_sample_info$normalized_name <- pc3_rankx2$normalized_name
pc3_sample_info$cell <- pc3_rankx2$cell

###########

RankedMatrix <- rbind(phh_rankx2, hepg2_rankx2)
RankedMatrix <- rbind(RankedMatrix, ha1e_rankx2)
RankedMatrix <- rbind(RankedMatrix, a375_rankx2)
RankedMatrix <- rbind(RankedMatrix, mcf7_rankx2)
RankedMatrix <- rbind(RankedMatrix, pc3_rankx2)

write.table(RankedMatrix, file="../output_files/RankedMatrix_x2.txt", quote=FALSE, sep="\t")

SampleInfoFile <- rbind(phh_sample_info, hepg2_sample_info)
SampleInfoFile <- rbind(SampleInfoFile, ha1e_sample_info)
SampleInfoFile <- rbind(SampleInfoFile, a375_sample_info)
SampleInfoFile <- rbind(SampleInfoFile, mcf7_sample_info)
SampleInfoFile <- rbind(SampleInfoFile, pc3_sample_info)

write.table(SampleInfoFile, file="../output_files/SampleInfoFile_x2.txt", quote=FALSE, sep="\t")


#########################################
# Kruskal Algorithm - Borda Merging 
#Kru-Bor merge samples by drug
#########################################
#RankedMatrix = fread.Ranked.Matrix(MatrixPath)
#SampleInfoFile = fread(SampleInfoFile, 
#                       data.table = FALSE)
#DrugList = read.Drugs(DrugFile)
#DrugEset = drug.eset(RankedMatrix, SampleInfoFile, DrugList)
#
#mergedset = RankMerging(DrugEset)
SampleInfoFile <- read.table("../output_files/SampleInfoFile_250.txt", sep="\t")

drug.eset<- function(RankedMatrix, SampleInfoFile, DrugList){
    #Identify samples treated with drugs of interest
    DrugInstances = SampleInfoFile$instance_id[which(SampleInfoFile$normalized_name%in%DrugList)]
    #make a subset of the RankedMatrix, containing only these Ids
    Drug.RankMatrix = RankedMatrix[,as.character(DrugInstances)]
    #dummy eset, to force sample order
    eset = ExpressionSet(assayData = as.matrix(RankedMatrix))
    #make drug annotation object for the ExpressionSet
    DrugData = SampleInfoFile[which(SampleInfoFile$instance_id%in%colnames(exprs(eset))), c(1, 2)] #id, drugname
    rownames(DrugData) = DrugData[,1]
    DrugData = DrugData[,-1, drop = FALSE]
    DrugAFD = new("AnnotatedDataFrame", data = DrugData)
    #make the ExpressionSet
    eset = ExpressionSet(assayData = as.matrix(RankedMatrix), phenoData = DrugAFD)
    return(eset)
}

#intersect(SampleInfoFile$normalized_name,DrugList)
#test <- unique(SampleInfoFile$normalized_name)

#DrugList <- as.data.frame(targets$CAM_ID)
#colnames(DrugList) <- "normalized_name"


#Kai
DrugList <- RankedMatrix$normalized_name
DrugList <- unique(DrugList)
#colnames(SampleInfoFile)[2] <- "normalized_name"
#Ranked250Matrix_train <- Ranked250Matrix[which(Ranked250Matrix$normalized_name%in%Training_CMID$CAM_ID),]
#RankedMatrix<-t(Ranked250Matrix_train)
#RankedMatrix<-RankedMatrix[1:12328,]
#write.csv(RankedMatrix2,file="DrugSig250.csv",quote=F)
#RankedMatrix250<-read.csv("DrugSig250.csv",row.names=1)

#RankedMatrix <- RankedMatrix250[which(RankedMatrix250$)]

RankedMatrix <- RankedMatrix[,1:12328]

DrugEset2 = drug.eset(t(RankedMatrix), SampleInfoFile, DrugList)
mergedset = RankMerging(DrugEset2)
dim(exprs(mergedset))

#RankedMatrixTraining <- RankedMatrix[which(RankedMatrix$normalized_name%in%Training_CMID$CAM_ID),]
#DrugsListT <- RankedMatrixTraining$normalized_name
#colnames(DrugListT) <- "normalized_name"

#DrugEsetT = drug.eset(t(RankedMatrixTraining), SampleInfoFile, DrugListT)
#mergedTset = RankMerging(DrugEsetT)

#merged100set <- as.data.frame(t(exprs(mergedset)))
#merged100set[merged100set<=100] <- 1
#merged100set[merged100set>=12228] <- -1
#merged100set[merged100set>1] <- 0

##Testing if this removes columns incorrectly

merged100set_test <- as.data.frame(t(exprs(mergedset)))
merged100set_test[merged100set_test<=100] <- 1
merged100set_test[merged100set_test>=12228] <- 1
merged100set_test[merged100set_test!=1] <- 0

merged100setT_test <- merged100set_test[which(rownames(merged100set_test)%in%Training_CMID$CAM_ID),]
merged100setT0_test <- merged100setT_test[,-(which(colSums(merged100setT_test)==0))]
merged100setT0g <- merged100setT_test[,apply(merged100setT_test,2,function(x) !all(x==0))]

#write.table(merged100set, "Merged100CAM.txt", quote=F, sep="\t")
#merged250set <- as.data.frame(t(exprs(mergedset)))
#merged250set[merged250set<=250] <- 1
#merged250set[merged250set>=12078] <- -1
#merged250set[merged250set>1] <- 0

merged250set_test <- as.data.frame(t(exprs(mergedset)))
merged250set_test[merged250set_test<=250] <- 1
merged250set_test[merged250set_test>=12078] <- 1
merged250set_test[merged250set_test!=1] <- 0

merged250setT_test <- merged250set_test[which(rownames(merged250set_test)%in%Training_CMID$CAM_ID),]
merged250setT0_test <- merged250setT_test[,-(which(colSums(merged250setT_test)==0))]
merged250setT0g <- merged250setT_test[,apply(merged250setT_test,2,function(x) !all(x==0))]

#write.table(merged250set, "Merged250CAM.txt", quote=F, sep="\t")
#merged500set <- as.data.frame(t(exprs(mergedset)))
#merged500set[merged500set<=500] <- 1
#merged500set[merged500set>=11828] <- -1
#merged500set[merged500set>1] <- 0

merged500set_test <- as.data.frame(t(exprs(mergedset)))
merged500set_test[merged500set_test<=500] <- 1
merged500set_test[merged500set_test>=11828] <- 1
merged500set_test[merged500set_test!=1] <- 0

merged500setT_test <- merged500set_test[which(rownames(merged500set_test)%in%Training_CMID$CAM_ID),]
merged500setT0_test <- merged500setT_test[,-(which(colSums(merged500setT_test)==0))]
merged500setT0g <- merged500setT_test[,apply(merged500setT_test,2,function(x) !all(x==0))]

#write.table(merged500set, "Merged500CAM.txt", quote=F, sep="\t")
#merged1000set <- as.data.frame(t(exprs(mergedset)))
#merged1000set[merged1000set<=1000] <- 1
#merged1000set[merged1000set>=11328] <- -1
#merged1000set[merged1000set>1] <- 0

merged1000set_test <- as.data.frame(t(exprs(mergedset)))
merged1000set_test[merged1000set_test<=1000] <- 1
merged1000set_test[merged1000set_test>=11328] <- 1
merged1000set_test[merged1000set_test!=1] <- 0

merged1000setT_test <- merged1000set_test[which(rownames(merged1000set_test)%in%Training_CMID$CAM_ID),]
merged1000setT0_test <- merged1000setT_test[,-(which(colSums(merged1000setT_test)==0))]
merged1000setT0g <- merged1000setT_test[,apply(merged1000setT_test,2,function(x) !all(x==0))]

#write.table(merged1000set, "Merged1000CAM.txt", quote=F, sep="\t")

Training2 <- Training_CMID[which(Training_CMID$CAM_ID %in% mergedset$normalized_name),]

#merged100set0 <- merged100set[,-(which(colSums(merged100set)==0))]
merged100setT <- merged100set[which(rownames(merged100set)%in%Training_CMID$CAM_ID),]
#merged100setT0 <- merged100setT[,-(which(colSums(merged100setT)==0))]
merged100setT0 <- merged100setT[,apply(merged100setT,2,function(x) !all(x==0))]


#merged250setT0 <- merged250setT[,-(which(colSums(merged250setT)==0))]
merged250setT0 <- merged250setT[,apply(merged250setT,2,function(x) !all(x==0))]

merged500setT <- merged500set[which(rownames(merged500set)%in%Training_CMID$CAM_ID),]
#merged500setT0 <- merged500setT[,-(which(colSums(merged500setT)==0))]
merged500setT0 <- merged500setT[,apply(merged500setT,2,function(x) !all(x==0))]

merged1000setT <- merged1000set[which(rownames(merged1000set)%in%Training_CMID$CAM_ID),]
#merged1000setT0 <- merged1000setT[,-(which(colSums(merged1000setT)==0))]
merged1000setT0 <- merged1000setT[,apply(merged1000setT,2,function(x) !all(x==0))]

#mergedallset <- t(exprs(mergedset))
#mergedallsetT <- mergedallset[which(rownames(mergedallset)%in%Training_CMID$CAM_ID),]


test1 <- as.data.frame(rowSums(merged100setT))
test2 <- as.data.frame(colSums(merged100setT))
colnames(test2) <- "count"
test20 <- as.data.frame(test2[which(test2$count==0),])

DrugListZ <- ZMatrix$normalized_name
DrugListZ <- unique(DrugListZ)

#DrugEsetZ <- drug.eset(t(ZMatrix[,1:12328]), SampleInfoZ, DrugListZ)
DrugEsetZ <- drug.eset(t(ZMatrix_train), SampleInfoZ, DrugListZ)
mergedZset <- RankMerging(DrugEsetZ)



############
#100
dili100_1 <- as.data.frame(Training_CMID[,2])
#dili100_1 <- cbind(Training2$DILI1, merged100setT0)
#dili100_1 <- cbind(Training2$DILI1, merged100setT)
merged100setT0x <- as.data.frame(merged100setT0g)
merged100setT0x$CAM_ID <- rownames(merged100setT0g)
dili100_1 <- merge(Training2[,c(1,2)], merged100setT0x, by="CAM_ID")
rownames(dili100_1) <- dili100_1$CAM_ID

dili100_1 <- dili100_1[,c(-1)]
colnames(dili100_1)[1] <- "dili_1"
dili100_1$dili_1 <- paste0('x',dili100_1$dili_1)
dili100_1$dili_1 <- as.factor(dili100_1$dili_1)
dili100_1 <- as.data.frame(dili100_1)

targ1 <- dili100_1$dili_1

dili100_1fish <- apply(dili100_1[,c(2:ncol(dili100_1))], 2, function (x) 
    fisher.test(table(targ1,x))$p.value)

dili100_1fish <- as.data.frame(dili100_1fish)
dili100_1fish$Gene <- rownames(dili100_1fish)
GeneList <- as.data.frame(dili100_1fish[which(dili100_1fish$dili100_1fish<=0.01),])
dili100_1g <- dili100_1[,names(dili100_1[,1:ncol(dili100_1)])%in%c(GeneList$Gene, "dili_1")]

GeneListx <- as.data.frame(dili100_1fish[which(dili100_1fish$dili100_1fish<0.05),])
dili100_1gx <- dili100_1[,names(dili100_1[,1:ncol(dili100_1)])%in%c(GeneListx$Gene, "dili_1")]


#set.seed(123)
#down100_1gx <- downSample(dili100_1gx[,!names(dili100_1gx)%in%c("dili_1")], as.factor(dili100_1$dili_1), yname ="dili_1")


#merged100setT0x <- as.data.frame(merged100setT0g)
#merged100setT0x$CAM_ID <- rownames(merged100setT0g)

#dili100_3 <- as.data.frame(Training_CMID[,4])
#dili100_3 <- cbind(Training2$DILI3, merged100setT0)
dili100_3 <- merge(Training2[,c(1,4)], merged100setT0x, by="CAM_ID")
rownames(dili100_3) <- dili100_3$CAM_ID
dili100_3 <- dili100_3[,c(-1)]
colnames(dili100_3)[1] <- "dili_3"
dili100_3 <- as.data.frame(dili100_3)
dili100_3$dili_3 <- paste0('x',dili100_3$dili_3)
dili100_3$dili_3 <- as.factor(dili100_3$dili_3)

targ3 <- dili100_3$dili_3

dili100_3fish <- apply(dili100_3[,c(2:ncol(dili100_3))], 2, function (x) 
    fisher.test(table(targ3,x))$p.value)

#dili100_1chi <- apply(dili100_1[,c(2:ncol(dili100_1))], 2, function (x) 
#  chisq.test(table(targ,x))$p.value)

dili100_3fish <- as.data.frame(dili100_3fish)
dili100_3fish$Gene <- rownames(dili100_3fish)
GeneList <- as.data.frame(dili100_3fish[which(dili100_3fish$dili100_3fish<=0.01),])
dili100_3g <- dili100_3[,names(dili100_3[,1:ncol(dili100_3)])%in%c(GeneList$Gene, "dili_3")]

#merged100setT0x <- as.data.frame(merged100setT0g)
#merged100setT0x$CAM_ID <- rownames(merged100setT0g)

#dili100_5 <- as.data.frame(Training_CMID[,6])
#dili100_5 <- cbind(Training2$DILI3, merged100setT0)
dili100_5 <- merge(Training2[,c(1,6)], merged100setT0x, by="CAM_ID")
rownames(dili100_5) <- dili100_5$CAM_ID
dili100_5 <- dili100_5[,c(-1)]
colnames(dili100_5)[1] <- "dili_5"
dili100_5 <- as.data.frame(dili100_5)
dili100_5$dili_5 <- paste0('x',dili100_5$dili_5)
dili100_5$dili_5 <- as.factor(dili100_5$dili_5)

targ5 <- dili100_5$dili_5

dili100_5fish <- apply(dili100_5[,c(2:ncol(dili100_5))], 2, function (x) 
    fisher.test(table(targ5,x))$p.value)

#dili100_1chi <- apply(dili100_1[,c(2:ncol(dili100_1))], 2, function (x) 
#  chisq.test(table(targ,x))$p.value)

dili100_5fish <- as.data.frame(dili100_5fish)
dili100_5fish$Gene <- rownames(dili100_5fish)
GeneList <- as.data.frame(dili100_5fish[which(dili100_5fish$dili100_5fish<=0.01),])
dili100_5g <- dili100_5[,names(dili100_5[,1:ncol(dili100_5)])%in%c(GeneList$Gene, "dili_5")]


#
#merged100setT0x <- as.data.frame(merged100setT0g)
#merged100setT0x$CAM_ID <- rownames(merged100setT0g)

#dili100_6 <- as.data.frame(Training_CMID[,6])
#dili100_6 <- cbind(Training2$DILI3, merged100setT0)
dili100_6 <- merge(Training2[,c(1,7)], merged100setT0x, by="CAM_ID")
rownames(dili100_6) <- dili100_6$CAM_ID
dili100_6 <- dili100_6[,c(-1)]
colnames(dili100_6)[1] <- "dili_6"
dili100_6 <- as.data.frame(dili100_6)
dili100_6$dili_6 <- paste0('x',dili100_6$dili_6)
dili100_6$dili_6 <- as.factor(dili100_6$dili_6)

targ6 <- dili100_6$dili_6

dili100_6fish <- apply(dili100_6[,c(2:ncol(dili100_6))], 2, function (x) 
    fisher.test(table(targ6,x))$p.value)

#dili100_1chi <- apply(dili100_1[,c(2:ncol(dili100_1))], 2, function (x) 
#  chisq.test(table(targ,x))$p.value)

dili100_6fish <- as.data.frame(dili100_6fish)
dili100_6fish$Gene <- rownames(dili100_6fish)
GeneList <- as.data.frame(dili100_6fish[which(dili100_6fish$dili100_6fish<=0.01),])
dili100_6g <- dili100_6[,names(dili100_6[,1:ncol(dili100_6)])%in%c(GeneList$Gene, "dili_6")]


########
#250
merged250setT0x <- as.data.frame(merged250setT0g)
merged250setT0x$CAM_ID <- rownames(merged250setT0g)

dili250_1 <- as.data.frame(Training_CMID[,2])
dili250_1 <- merge(Training2[,c(1,2)], merged250setT0x, by="CAM_ID")
rownames(dili250_1) <- dili250_1$CAM_ID
dili250_1 <- dili250_1[,c(-1)]
colnames(dili250_1)[1] <- "dili_1"
dili250_1 <- as.data.frame(dili250_1)
dili250_1$dili_1 <- paste0('x',dili250_1$dili_1)
dili250_1$dili_1 <- as.factor(dili250_1$dili_1)

targ1 <- dili250_1$dili_1

dili250_1fish <- apply(dili250_1[,c(2:ncol(dili250_1))], 2, function (x) 
    fisher.test(table(targ1,x))$p.value)

#dili100_1chi <- apply(dili100_1[,c(2:ncol(dili100_1))], 2, function (x) 
#  chisq.test(table(targ,x))$p.value)

dili250_1fish <- as.data.frame(dili250_1fish)
dili250_1fish$Gene <- rownames(dili250_1fish)
GeneList <- as.data.frame(dili250_1fish[which(dili250_1fish$dili250_1fish<=0.01),])
dili250_1g <- dili250_1[,names(dili250_1[,1:ncol(dili250_1)])%in%c(GeneList$Gene, "dili_1")]


dili250_3 <- merge(Training2[,c(1,4)], merged250setT0x, by="CAM_ID")
rownames(dili250_3) <- dili250_3$CAM_ID
dili250_3 <- dili250_3[,c(-1)]
colnames(dili250_3)[1] <- "dili_3"
dili250_3 <- as.data.frame(dili250_3)
dili250_3$dili_3 <- paste0('x',dili250_3$dili_3)
dili250_3$dili_3 <- as.factor(dili250_3$dili_3)

targ3 <- dili250_3$dili_3
dili250_3fish <- apply(dili250_3[,c(2:ncol(dili250_3))], 2, function (x) 
    fisher.test(table(targ3,x))$p.value)

dili250_3fish <- as.data.frame(dili250_3fish)
dili250_3fish$Gene <- rownames(dili250_3fish)
GeneList <- as.data.frame(dili250_3fish[which(dili250_3fish$dili250_3fish<=0.01),])
dili250_3g <- dili250_3[,names(dili250_3[,1:ncol(dili250_3)])%in%c(GeneList$Gene, "dili_3")]


dili250_5 <- as.data.frame(Training_CMID[,2])
dili250_5 <- merge(Training2[,c(1,6)], merged250setT0x, by="CAM_ID")
rownames(dili250_5) <- dili250_5$CAM_ID
dili250_5 <- dili250_5[,c(-1)]
colnames(dili250_5)[1] <- "dili_5"
dili250_5 <- as.data.frame(dili250_5)
dili250_5$dili_5 <- paste0('x',dili250_5$dili_5)
dili250_5$dili_5 <- as.factor(dili250_5$dili_5)

targ5 <- dili250_5$dili_5
dili250_5fish <- apply(dili250_5[,c(2:ncol(dili250_5))], 2, function (x) 
    fisher.test(table(targ5,x))$p.value)

dili250_5fish <- as.data.frame(dili250_5fish)
dili250_5fish$Gene <- rownames(dili250_5fish)
GeneList <- as.data.frame(dili250_5fish[which(dili250_5fish$dili250_5fish<=0.01),])
dili250_5g <- dili250_5[,names(dili250_5[,1:ncol(dili250_5)])%in%c(GeneList$Gene, "dili_5")]

dili250_6 <- as.data.frame(Training_CMID[,2])
dili250_6 <- merge(Training2[,c(1,7)], merged250setT0x, by="CAM_ID")
rownames(dili250_6) <- dili250_6$CAM_ID
dili250_6 <- dili250_6[,c(-1)]
colnames(dili250_6)[1] <- "dili_6"
dili250_6 <- as.data.frame(dili250_6)
dili250_6$dili_6 <- paste0('x',dili250_6$dili_6)
dili250_6$dili_6 <- as.factor(dili250_6$dili_6)

targ6 <- dili250_6$dili_6
dili250_6fish <- apply(dili250_6[,c(2:ncol(dili250_6))], 2, function (x) 
    fisher.test(table(targ6,x))$p.value)

dili250_6fish <- as.data.frame(dili250_6fish)
dili250_6fish$Gene <- rownames(dili250_6fish)
GeneList <- as.data.frame(dili250_6fish[which(dili250_6fish$dili250_6fish<=0.01),])
dili250_6g <- dili250_6[,names(dili250_6[,1:ncol(dili250_6)])%in%c(GeneList$Gene, "dili_6")]


#######
merged500setT0x <- as.data.frame(merged500setT0g)
merged500setT0x$CAM_ID <- rownames(merged500setT0g)

##500
dili500_1 <- merge(Training2[,c(1,2)], merged500setT0x, by="CAM_ID")
rownames(dili500_1) <- dili500_1$CAM_ID
dili500_1 <- dili500_1[,c(-1)]
colnames(dili500_1)[1] <- "dili_1"
dili500_1 <- as.data.frame(dili500_1)
dili500_1$dili_1 <- paste0('x',dili500_1$dili_1)
dili500_1$dili_1 <- as.factor(dili500_1$dili_1)

targ1 <- dili500_1$dili_1

dili500_1fish <- apply(dili500_1[,c(2:ncol(dili500_1))], 2, function (x) 
    fisher.test(table(targ1,x))$p.value)

#dili100_1chi <- apply(dili100_1[,c(2:ncol(dili100_1))], 2, function (x) 
#  chisq.test(table(targ,x))$p.value)

dili500_1fish <- as.data.frame(dili500_1fish)
dili500_1fish$Gene <- rownames(dili500_1fish)
GeneList <- as.data.frame(dili500_1fish[which(dili500_1fish$dili500_1fish<=0.01),])
dili500_1g <- dili500_1[,names(dili500_1[,1:ncol(dili500_1)])%in%c(GeneList$Gene, "dili_1")]


###
test <- merge(dili250_1g, dili500_1g, by="")

###

dili500_3 <- merge(Training2[,c(1,4)], merged500setT0x, by="CAM_ID")
rownames(dili500_3) <- dili500_3$CAM_ID
dili500_3 <- dili500_3[,c(-1)]
colnames(dili500_3)[1] <- "dili_3"
dili500_3 <- as.data.frame(dili500_3)
dili500_3$dili_3 <- paste0('x',dili500_3$dili_3)
dili500_3$dili_3 <- as.factor(dili500_3$dili_3)

targ3 <- dili500_3$dili_3
dili500_3fish <- apply(dili500_3[,c(2:ncol(dili500_3))], 2, function (x) 
    fisher.test(table(targ3,x))$p.value)

dili500_3fish <- as.data.frame(dili500_3fish)
dili500_3fish$Gene <- rownames(dili500_3fish)
GeneList <- as.data.frame(dili500_3fish[which(dili500_3fish$dili500_3fish<=0.01),])
dili500_3g <- dili500_3[,names(dili500_3[,1:ncol(dili500_3)])%in%c(GeneList$Gene, "dili_3")]

dili500_5 <- merge(Training2[,c(1,6)], merged500setT0x, by="CAM_ID")
rownames(dili500_5) <- dili500_5$CAM_ID
dili500_5 <- dili500_5[,c(-1)]
colnames(dili500_5)[1] <- "dili_5"
dili500_5 <- as.data.frame(dili500_5)
dili500_5$dili_5 <- paste0('x',dili500_5$dili_5)
dili500_5$dili_5 <- as.factor(dili500_5$dili_5)

targ5 <- dili500_5$dili_5
dili500_5fish <- apply(dili500_5[,c(2:ncol(dili500_5))], 2, function (x) 
    fisher.test(table(targ5,x))$p.value)

dili500_5fish <- as.data.frame(dili500_5fish)
dili500_5fish$Gene <- rownames(dili500_5fish)
GeneList <- as.data.frame(dili500_5fish[which(dili500_5fish$dili500_5fish<=0.01),])
dili500_5g <- dili500_5[,names(dili500_5[,1:ncol(dili500_5)])%in%c(GeneList$Gene, "dili_5")]

dili500_6 <- merge(Training2[,c(1,7)], merged500setT0x, by="CAM_ID")
rownames(dili500_6) <- dili500_6$CAM_ID
dili500_6 <- dili500_6[,c(-1)]
colnames(dili500_6)[1] <- "dili_6"
dili500_6 <- as.data.frame(dili500_6)
dili500_6$dili_6 <- paste0('x',dili500_6$dili_6)
dili500_6$dili_6 <- as.factor(dili500_6$dili_6)

targ6 <- dili500_6$dili_6
dili500_6fish <- apply(dili500_6[,c(2:ncol(dili500_6))], 2, function (x) 
    fisher.test(table(targ6,x))$p.value)

dili500_6fish <- as.data.frame(dili500_6fish)
dili500_6fish$Gene <- rownames(dili500_6fish)
GeneList <- as.data.frame(dili500_6fish[which(dili500_6fish$dili500_6fish<=0.01),])
dili500_6g <- dili500_6[,names(dili500_6[,1:ncol(dili500_6)])%in%c(GeneList$Gene, "dili_6")]
######
merged1000setT0x <- as.data.frame(merged1000setT0g)
merged1000setT0x$CAM_ID <- rownames(merged1000setT0g)

##1000
dili1000_1 <- merge(Training2[,c(1,2)], merged1000setT0x, by="CAM_ID")
rownames(dili1000_1) <- dili1000_1$CAM_ID
dili1000_1 <- dili1000_1[,c(-1)]
colnames(dili1000_1)[1] <- "dili_1"
dili1000_1 <- as.data.frame(dili1000_1)
dili1000_1$dili_1 <- paste0('x',dili1000_1$dili_1)
dili1000_1$dili_1 <- as.factor(dili1000_1$dili_1)

targ1 <- dili1000_1$dili_1

dili1000_1fish <- apply(dili1000_1[,c(2:ncol(dili1000_1))], 2, function (x) 
    fisher.test(table(targ1,x))$p.value)

#dili100_1chi <- apply(dili100_1[,c(2:ncol(dili100_1))], 2, function (x) 
#  chisq.test(table(targ,x))$p.value)

dili1000_1fish <- as.data.frame(dili1000_1fish)
dili1000_1fish$Gene <- rownames(dili1000_1fish)
GeneList <- as.data.frame(dili1000_1fish[which(dili1000_1fish$dili1000_1fish<=0.01),])
dili1000_1g <- dili1000_1[,names(dili1000_1[,1:ncol(dili1000_1)])%in%c(GeneList$Gene, "dili_1")]

dili1000_3 <- merge(Training2[,c(1,4)], merged1000setT0x, by="CAM_ID")
rownames(dili1000_3) <- dili1000_3$CAM_ID
dili1000_3 <- dili1000_3[,c(-1)]
colnames(dili1000_3)[1] <- "dili_3"
dili1000_3 <- as.data.frame(dili1000_3)
dili1000_3$dili_3 <- paste0('x',dili1000_3$dili_3)
dili1000_3$dili_3 <- as.factor(dili1000_3$dili_3)

targ3 <- dili1000_3$dili_3
dili1000_3fish <- apply(dili1000_3[,c(2:ncol(dili1000_3))], 2, function (x) 
    fisher.test(table(targ3,x))$p.value)

dili1000_3fish <- as.data.frame(dili1000_3fish)
dili1000_3fish$Gene <- rownames(dili1000_3fish)
GeneList <- as.data.frame(dili1000_3fish[which(dili1000_3fish$dili1000_3fish<=0.01),])
dili1000_3g <- dili1000_3[,names(dili1000_3[,1:ncol(dili1000_3)])%in%c(GeneList$Gene, "dili_3")]

dili1000_5 <- merge(Training2[,c(1,6)], merged1000setT0x, by="CAM_ID")
rownames(dili1000_5) <- dili1000_5$CAM_ID
dili1000_5 <- dili1000_5[,c(-1)]
colnames(dili1000_5)[1] <- "dili_5"
dili1000_5 <- as.data.frame(dili1000_5)
dili1000_5$dili_5 <- paste0('x',dili1000_5$dili_5)
dili1000_5$dili_5 <- as.factor(dili1000_5$dili_5)

targ5 <- dili1000_5$dili_5
dili1000_5fish <- apply(dili1000_5[,c(2:ncol(dili1000_5))], 2, function (x) 
    fisher.test(table(targ5,x))$p.value)

dili1000_5fish <- as.data.frame(dili1000_5fish)
dili1000_5fish$Gene <- rownames(dili1000_5fish)
GeneList <- as.data.frame(dili1000_5fish[which(dili1000_5fish$dili1000_5fish<=0.01),])
dili1000_5g <- dili1000_5[,names(dili1000_5[,1:ncol(dili1000_5)])%in%c(GeneList$Gene, "dili_5")]

dili1000_6 <- merge(Training2[,c(1,7)], merged1000setT0x, by="CAM_ID")
rownames(dili1000_6) <- dili1000_6$CAM_ID
dili1000_6 <- dili1000_6[,c(-1)]
colnames(dili1000_6)[1] <- "dili_6"
dili1000_6 <- as.data.frame(dili1000_6)
dili1000_6$dili_6 <- paste0('x',dili1000_6$dili_6)
dili1000_6$dili_6 <- as.factor(dili1000_6$dili_6)

targ6 <- dili1000_6$dili_6
dili1000_6fish <- apply(dili1000_6[,c(2:ncol(dili1000_6))], 2, function (x) 
    fisher.test(table(targ6,x))$p.value)

dili1000_6fish <- as.data.frame(dili1000_6fish)
dili1000_6fish$Gene <- rownames(dili1000_6fish)
GeneList <- as.data.frame(dili1000_6fish[which(dili1000_6fish$dili1000_6fish<=0.01),])
dili1000_6g <- dili1000_6[,names(dili1000_6[,1:ncol(dili1000_6)])%in%c(GeneList$Gene, "dili_6")]

############
cl <- makeCluster(30)

registerDoParallel(cl)

####
ctl <- trainControl(method="repeatedcv", number=5,repeats = 100, returnResamp="all",allowParallel = T,savePredictions = T,
                    classProbs=TRUE, summaryFunction=twoClassSummary)

##Testing Signature Sizing
set.seed(123)
rf100.3g <- train(dili_3~.,data=dili100_3g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf100.5g <- train(dili_5~.,data=dili100_5g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)

set.seed(123)
#11k genes
#rf100.1g <- train(dili_1~.,data=dili100_1, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
#0.01
rf100.1g <- train(dili_1~.,data=dili100_1g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
#0.05
#rf100.1gx <- train(dili_1~.,data=dili100_1gx, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
#downsampled
#rf100down.1gx <- train(dili_1~.,data=down100_1gx, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
#0.01
rf100.6g <- train(dili_6~.,data=dili100_6g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)

set.seed(123)
svm100.1 <- train(dili_1~., data=dili100_1g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm100.3 <- train(dili_3~., data=dili100_3g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm100.5 <- train(dili_5~., data=dili100_5g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm100.6 <- train(dili_6~., data=dili100_6g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
#
set.seed(123)
rpart100.1 <- train(dili_1~., data=dili100_1g, method="rpart", trControl=ctl, metric = "ROC")
rpart100.3 <- train(dili_3~., data=dili100_3g, method="rpart", trControl=ctl, metric = "ROC")
rpart100.5 <- train(dili_5~., data=dili100_5g, method="rpart", trControl=ctl, metric = "ROC")
rpart100.6 <- train(dili_6~., data=dili100_6g, method="rpart", trControl=ctl, metric = "ROC")

set.seed(123)
glm100.1 <- train(dili_1~., data=dili100_1g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm100.3 <- train(dili_3~., data=dili100_3g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm100.5 <- train(dili_5~., data=dili100_5g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm100.6 <- train(dili_6~., data=dili100_6g, method="glm", trControl = ctl, metric = "ROC", family="binomial")

set.seed(123)
NB100.1 <- train(dili_1~., data=dili100_1g, method="naive_bayes", trControl = ctl, metric="ROC")
NB100.3 <- train(dili_3~., data=dili100_3g, method="naive_bayes", trControl = ctl, metric="ROC")
NB100.5 <- train(dili_5~., data=dili100_5g, method="naive_bayes", trControl = ctl, metric="ROC")
NB100.6 <- train(dili_6~., data=dili100_6g, method="naive_bayes", trControl = ctl, metric="ROC")

#stopCluster(cl)


#
set.seed(123)
rf250.1 <- train(dili_1~.,data=dili250_1g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf250.3 <- train(dili_3~.,data=dili250_3g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf250.5 <- train(dili_5~.,data=dili250_5g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf250.6 <- train(dili_6~.,data=dili250_6g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)

set.seed(123)
svm250.1 <- train(dili_1~., data=dili250_1g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm250.3 <- train(dili_3~., data=dili250_3g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm250.5 <- train(dili_5~., data=dili250_5g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm250.6 <- train(dili_6~., data=dili250_6g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)

set.seed(123)
rpart250.1 <- train(dili_1~., data=dili250_1g, method="rpart", trControl=ctl, metric = "ROC")
rpart250.3 <- train(dili_3~., data=dili250_3g, method="rpart", trControl=ctl, metric = "ROC")
rpart250.5 <- train(dili_5~., data=dili250_5g, method="rpart", trControl=ctl, metric = "ROC")
rpart250.6 <- train(dili_6~., data=dili250_6g, method="rpart", trControl=ctl, metric = "ROC")

set.seed(123)
glm250.1 <- train(dili_1~., data=dili250_1g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm250.3 <- train(dili_3~., data=dili250_3g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm250.5 <- train(dili_5~., data=dili250_5g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm250.6 <- train(dili_6~., data=dili250_6g, method="glm", trControl = ctl, metric = "ROC", family="binomial")

set.seed(123)
NB250.1 <- train(dili_1~., data=dili250_1g, method="naive_bayes", trControl = ctl, metric="ROC")
NB250.3 <- train(dili_3~., data=dili250_3g, method="naive_bayes", trControl = ctl, metric="ROC")
NB250.5 <- train(dili_5~., data=dili250_5g, method="naive_bayes", trControl = ctl, metric="ROC")
NB250.6 <- train(dili_6~., data=dili250_6g, method="naive_bayes", trControl = ctl, metric="ROC")
#
set.seed(123)
rf500.1 <- train(dili_1~.,data=dili500_1g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf500.3 <- train(dili_3~.,data=dili500_3g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf500.5 <- train(dili_5~.,data=dili500_5g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf500.6 <- train(dili_6~.,data=dili500_6g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)

set.seed(123)
svm500.1 <- train(dili_1~., data=dili500_1g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm500.3 <- train(dili_3~., data=dili500_3g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm500.5 <- train(dili_5~., data=dili500_5g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm500.6 <- train(dili_6~., data=dili500_6g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)

set.seed(123)
rpart500.1 <- train(dili_1~., data=dili500_1g, method="rpart", trControl=ctl, metric = "ROC")
rpart500.3 <- train(dili_3~., data=dili500_3g, method="rpart", trControl=ctl, metric = "ROC")
rpart500.5 <- train(dili_5~., data=dili500_5g, method="rpart", trControl=ctl, metric = "ROC")
rpart500.6 <- train(dili_6~., data=dili500_6g, method="rpart", trControl=ctl, metric = "ROC")

set.seed(123)
glm500.1 <- train(dili_1~., data=dili500_1g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm500.3 <- train(dili_3~., data=dili500_3g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm500.5 <- train(dili_5~., data=dili500_5g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm500.6 <- train(dili_6~., data=dili500_6g, method="glm", trControl = ctl, metric = "ROC", family="binomial")

set.seed(123)
NB500.1 <- train(dili_1~., data=dili500_1g, method="naive_bayes", trControl = ctl, metric="ROC")
NB500.3 <- train(dili_3~., data=dili500_3g, method="naive_bayes", trControl = ctl, metric="ROC")
NB500.5 <- train(dili_5~., data=dili500_5g, method="naive_bayes", trControl = ctl, metric="ROC")
NB500.6 <- train(dili_6~., data=dili500_6g, method="naive_bayes", trControl = ctl, metric="ROC")
#
set.seed(123)
rf1000.1 <- train(dili_1~.,data=dili1000_1g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf1000.3 <- train(dili_3~.,data=dili1000_3g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf1000.5 <- train(dili_5~.,data=dili1000_5g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)
rf1000.6 <- train(dili_6~.,data=dili1000_6g, method="rf", trControl=ctl, metric="ROC", allowParallel=T)

set.seed(123)
svm1000.1 <- train(dili_1~., data=dili1000_1g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm1000.3 <- train(dili_3~., data=dili1000_3g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm1000.5 <- train(dili_5~., data=dili1000_5g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)
svm1000.6 <- train(dili_6~., data=dili1000_6g, method="svmRadial", trControl = ctl, metric ="ROC", allowParallel=T)

set.seed(123)
rpart1000.1 <- train(dili_1~., data=dili1000_1g, method="rpart", trControl=ctl, metric = "ROC")
rpart1000.3 <- train(dili_3~., data=dili1000_3g, method="rpart", trControl=ctl, metric = "ROC")
rpart1000.5 <- train(dili_5~., data=dili1000_5g, method="rpart", trControl=ctl, metric = "ROC")
rpart1000.6 <- train(dili_6~., data=dili1000_6g, method="rpart", trControl=ctl, metric = "ROC")

set.seed(123)
glm1000.1 <- train(dili_1~., data=dili1000_1g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm1000.3 <- train(dili_3~., data=dili1000_3g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm1000.5 <- train(dili_5~., data=dili1000_5g, method="glm", trControl = ctl, metric = "ROC", family="binomial")
glm1000.6 <- train(dili_6~., data=dili1000_6g, method="glm", trControl = ctl, metric = "ROC", family="binomial")

set.seed(123)
NB1000.1 <- train(dili_1~., data=dili1000_1g, method="naive_bayes", trControl = ctl, metric="ROC")
NB1000.3 <- train(dili_3~., data=dili1000_3g, method="naive_bayes", trControl = ctl, metric="ROC")
NB1000.5 <- train(dili_5~., data=dili1000_5g, method="naive_bayes", trControl = ctl, metric="ROC")
NB1000.6 <- train(dili_6~., data=dili1000_6g, method="naive_bayes", trControl = ctl, metric="ROC")

stopCluster(cl)

################
#TrainCMIDs for rf 100 merge
#####
trmerge100_dili1_svm <- predict(svm100.1, merged100setT)
names(trmerge100_dili1_svm) <- rownames(merged100setT)
trmerge100_dili1_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge100_dili1_svm2) <- "CMID"
trmerge100_dili1_svm2$Dili_1 <- as.character(trmerge100_dili1_svm[Training_CMID$CAM_ID])

trmerge100_dili1_svm2$Dili_1 <- sub(".", "", trmerge100_dili1_svm2$Dili_1)

write.csv(trmerge100_dili1_svm2, "../output_files/Dili1_TrainingPrediction_merge100.csv", quote=FALSE)
#
trmerge250_dili1_svm <- predict(svm250.1, merged250setT)
names(trmerge250_dili1_svm) <- rownames(merged250setT)
trmerge250_dili1_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge250_dili1_svm2) <- "CMID"
trmerge250_dili1_svm2$Dili_1 <- as.character(trmerge250_dili1_svm[Training_CMID$CAM_ID])

trmerge250_dili1_svm2$Dili_1 <- sub(".", "", trmerge250_dili1_svm2$Dili_1)

write.csv(trmerge250_dili1_svm2, "../output_files/Dili1_TrainingPrediction_merge250.csv", quote=FALSE)
#
trmerge500_dili1_svm <- predict(svm500.1, merged500setT)
names(trmerge500_dili1_svm) <- rownames(merged500setT)
trmerge500_dili1_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge500_dili1_svm2) <- "CMID"
trmerge500_dili1_svm2$Dili_1 <- as.character(trmerge500_dili1_svm[Training_CMID$CAM_ID])

trmerge500_dili1_svm2$Dili_1 <- sub(".", "", trmerge500_dili1_svm2$Dili_1)

write.csv(trmerge500_dili1_svm2, "../output_files/Dili1_TrainingPrediction_merge500.csv", quote=FALSE)
#
trmerge1000_dili1_svm <- predict(svm1000.1, merged1000setT)
names(trmerge1000_dili1_svm) <- rownames(merged1000setT)
trmerge1000_dili1_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge1000_dili1_svm2) <- "CMID"
trmerge1000_dili1_svm2$Dili_1 <- as.character(trmerge1000_dili1_svm[Training_CMID$CAM_ID])

trmerge1000_dili1_svm2$Dili_1 <- sub(".", "", trmerge1000_dili1_svm2$Dili_1)

write.csv(trmerge1000_dili1_svm2, "../output_files/Dili1_TrainingPrediction_merge1000.csv", quote=FALSE)
####

trmerge100_dili3_svm <- predict(svm100.3, merged100setT)
names(trmerge100_dili3_svm) <- rownames(merged100setT)
trmerge100_dili3_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge100_dili3_svm2) <- "CMID"
trmerge100_dili3_svm2$Dili_3 <- as.character(trmerge100_dili3_svm[Training_CMID$CAM_ID])

trmerge100_dili3_svm2$Dili_3 <- sub(".", "", trmerge100_dili3_svm2$Dili_3)

write.csv(trmerge100_dili3_svm2, "../output_files/Dili3_TrainingPrediction_merge100.csv", quote=FALSE)
#
trmerge250_dili3_svm <- predict(svm250.3, merged250setT)
names(trmerge250_dili3_svm) <- rownames(merged250setT)
trmerge250_dili3_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge250_dili3_svm2) <- "CMID"
trmerge250_dili3_svm2$Dili_3 <- as.character(trmerge250_dili3_svm[Training_CMID$CAM_ID])

trmerge250_dili3_svm2$Dili_3 <- sub(".", "", trmerge250_dili3_svm2$Dili_3)

write.csv(trmerge250_dili3_svm2, "../output_files/Dili3_TrainingPrediction_merge250.csv", quote=FALSE)
#
trmerge500_dili3_svm <- predict(svm500.3, merged500setT)
names(trmerge500_dili3_svm) <- rownames(merged500setT)
trmerge500_dili3_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge500_dili3_svm2) <- "CMID"
trmerge500_dili3_svm2$Dili_3 <- as.character(trmerge500_dili3_svm[Training_CMID$CAM_ID])

trmerge500_dili3_svm2$Dili_3 <- sub(".", "", trmerge500_dili3_svm2$Dili_3)

write.csv(trmerge500_dili3_svm2, "../output_files/Dili3_TrainingPrediction_merge500.csv", quote=FALSE)
#
trmerge1000_dili3_svm <- predict(svm1000.3, merged1000setT)
names(trmerge1000_dili3_svm) <- rownames(merged1000setT)
trmerge1000_dili3_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge1000_dili3_svm2) <- "CMID"
trmerge1000_dili3_svm2$Dili_3 <- as.character(trmerge1000_dili3_svm[Training_CMID$CAM_ID])

trmerge1000_dili3_svm2$Dili_3 <- sub(".", "", trmerge1000_dili3_svm2$Dili_3)

write.csv(trmerge1000_dili3_svm2, "../output_files/Dili3_TrainingPrediction_merge1000.csv", quote=FALSE)
####
trmerge100_dili5_svm <- predict(svm100.5, merged100setT)
names(trmerge100_dili5_svm) <- rownames(merged100setT)
trmerge100_dili5_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge100_dili5_svm2) <- "CMID"
trmerge100_dili5_svm2$Dili_5 <- as.character(trmerge100_dili5_svm[Training_CMID$CAM_ID])

trmerge100_dili5_svm2$Dili_5 <- sub(".", "", trmerge100_dili5_svm2$Dili_5)

write.csv(trmerge100_dili5_svm2, "../output_files/Dili5_TrainingPrediction_merge100.csv", quote=FALSE)
#
trmerge250_dili5_svm <- predict(svm250.1, merged250setT)
names(trmerge250_dili5_svm) <- rownames(merged250setT)
trmerge250_dili5_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge250_dili5_svm2) <- "CMID"
trmerge250_dili5_svm2$Dili_5 <- as.character(trmerge250_dili5_svm[Training_CMID$CAM_ID])

trmerge250_dili5_svm2$Dili_5 <- sub(".", "", trmerge250_dili5_svm2$Dili_5)

write.csv(trmerge250_dili5_svm2, "../output_files/Dili5_TrainingPrediction_merge250.csv", quote=FALSE)
#
trmerge500_dili5_svm <- predict(svm500.5, merged500setT)
names(trmerge500_dili5_svm) <- rownames(merged500setT)
trmerge500_dili5_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge500_dili5_svm2) <- "CMID"
trmerge500_dili5_svm2$Dili_5 <- as.character(trmerge500_dili5_svm[Training_CMID$CAM_ID])

trmerge500_dili5_svm2$Dili_5 <- sub(".", "", trmerge500_dili5_svm2$Dili_5)

write.csv(trmerge500_dili5_svm2, "../output_files/Dili5_TrainingPrediction_merge500.csv", quote=FALSE)
#
trmerge1000_dili5_svm <- predict(svm1000.5, merged1000setT)
names(trmerge1000_dili5_svm) <- rownames(merged1000setT)
trmerge1000_dili5_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge1000_dili5_svm2) <- "CMID"
trmerge1000_dili5_svm2$Dili_5 <- as.character(trmerge1000_dili5_svm[Training_CMID$CAM_ID])

trmerge1000_dili5_svm2$Dili_5 <- sub(".", "", trmerge1000_dili5_svm2$Dili_5)

write.csv(trmerge1000_dili5_svm2, "../output_files/Dili5_TrainingPrediction_merge1000.csv", quote=FALSE)
####
trmerge100_dili6_svm <- predict(svm100.6, merged100setT)
names(trmerge100_dili6_svm) <- rownames(merged100setT)
trmerge100_dili6_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge100_dili6_svm2) <- "CMID"
trmerge100_dili6_svm2$Dili_6 <- as.character(trmerge100_dili6_svm[Training_CMID$CAM_ID])

trmerge100_dili6_svm2$Dili_6 <- sub(".", "", trmerge100_dili6_svm2$Dili_6)

write.csv(trmerge100_dili6_svm2, "../output_files/Dili6_TrainingPrediction_merge100.csv", quote=FALSE)
#
trmerge250_dili6_svm <- predict(svm250.6, merged250setT)
names(trmerge250_dili6_svm) <- rownames(merged250setT)
trmerge250_dili6_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge250_dili6_svm2) <- "CMID"
trmerge250_dili6_svm2$Dili_6 <- as.character(trmerge250_dili6_svm[Training_CMID$CAM_ID])

trmerge250_dili6_svm2$Dili_6 <- sub(".", "", trmerge250_dili6_svm2$Dili_6)

write.csv(trmerge250_dili6_svm2, "../output_files/Dili6_TrainingPrediction_merge250.csv", quote=FALSE)
#
trmerge500_dili6_svm <- predict(svm500.6, merged500setT)
names(trmerge500_dili6_svm) <- rownames(merged500setT)
trmerge500_dili6_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge500_dili6_svm2) <- "CMID"
trmerge500_dili6_svm2$Dili_6 <- as.character(trmerge500_dili6_svm[Training_CMID$CAM_ID])

trmerge500_dili6_svm2$Dili_6 <- sub(".", "", trmerge500_dili6_svm2$Dili_6)

write.csv(trmerge500_dili6_svm2, "../output_files/Dili6_TrainingPrediction_merge500.csv", quote=FALSE)
#
trmerge1000_dili6_svm <- predict(svm1000.6, merged1000setT)
names(trmerge1000_dili6_svm) <- rownames(merged1000setT)
trmerge1000_dili6_svm2 <- as.data.frame(Training_CMID$CAM_ID, stringsAsFactors=FALSE)
colnames(trmerge1000_dili6_svm2) <- "CMID"
trmerge1000_dili6_svm2$Dili_6 <- as.character(trmerge1000_dili6_svm[Training_CMID$CAM_ID])

trmerge1000_dili6_svm2$Dili_6 <- sub(".", "", trmerge1000_dili6_svm2$Dili_6)

write.csv(trmerge1000_dili6_svm2, "../output_files/Dili6_TrainingPrediction_merge1000.csv", quote=FALSE)


###########
##Test Prediction Subset
#######
#TestCMIDs for rf 100 merge
Test <- targets[which(targets$Training_Validation==""),]

testmerge100_dili1_svm <- predict(svm100.1, merged100set)
names(testmerge100_dili1_svm) <- rownames(merged100set)
testmerge100_dili1_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili1_svm2) <- "CAM_ID"
testmerge100_dili1_svm2$Dili_1 <- as.character(testmerge100_dili1_svm[Test$CAM_ID])

testmerge100_dili1_svm2$Dili_1 <- sub(".", "", testmerge100_dili1_svm2$Dili_1)

#
testmerge250_dili1_svm <- predict(svm250.1, merged250set)
names(testmerge250_dili1_svm) <- rownames(merged250set)
testmerge250_dili1_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili1_svm2) <- "CAM_ID"
testmerge250_dili1_svm2$Dili_1 <- as.character(testmerge250_dili1_svm[Test$CAM_ID])

testmerge250_dili1_svm2$Dili_1 <- sub(".", "", testmerge250_dili1_svm2$Dili_1)

#write.csv(trmerge250_dili1_svm2, "Dili1_TrainingPrediction_merge250.csv", quote=FALSE)
#
testmerge500_dili1_svm <- predict(svm500.1, merged500set)
names(testmerge500_dili1_svm) <- rownames(merged500set)
testmerge500_dili1_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili1_svm2) <- "CAM_ID"
testmerge500_dili1_svm2$Dili_1 <- as.character(testmerge500_dili1_svm[Test$CAM_ID])

testmerge500_dili1_svm2$Dili_1 <- sub(".", "", testmerge500_dili1_svm2$Dili_1)

#write.csv(trmerge500_dili1_svm2, "Dili1_TrainingPrediction_merge500.csv", quote=FALSE)
#
testmerge1000_dili1_nb <- predict(NB1000.1, merged1000set)
names(testmerge1000_dili1_nb) <- rownames(merged1000set)
testmerge1000_dili1_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili1_nb2) <- "CAM_ID"
testmerge1000_dili1_nb2$Dili_1 <- as.character(testmerge1000_dili1_nb[Test$CAM_ID])

testmerge1000_dili1_nb2$Dili_1 <- sub(".", "", testmerge1000_dili1_nb2$Dili_1)

#write.csv(trmerge1000_dili1_svm2, "Dili1_TrainingPrediction_merge1000.csv", quote=FALSE)
####

testmerge100_dili3_rf <- predict(rf100.3g, merged100set)
names(testmerge100_dili3_rf) <- rownames(merged100set)
testmerge100_dili3_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili3_rf2) <- "CAM_ID"
testmerge100_dili3_rf2$Dili_3 <- as.character(testmerge100_dili3_rf[Test$CAM_ID])

testmerge100_dili3_rf2$Dili_3 <- sub(".", "", testmerge100_dili3_rf2$Dili_3)

#write.csv(trmerge100_dili3_svm2, "Dili3_TrainingPrediction_merge100.csv", quote=FALSE)
#
testmerge250_dili3_svm <- predict(svm250.3, merged250set)
names(testmerge250_dili3_svm) <- rownames(merged250set)
testmerge250_dili3_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili3_svm2) <- "CAM_ID"
testmerge250_dili3_svm2$Dili_3 <- as.character(testmerge250_dili3_svm[Test$CAM_ID])

testmerge250_dili3_svm2$Dili_3 <- sub(".", "", testmerge250_dili3_svm2$Dili_3)

#write.csv(trmerge250_dili3_svm2, "Dili3_TrainingPrediction_merge250.csv", quote=FALSE)
#
testmerge250_dili3_rf <- predict(rf250.3, merged250set)
names(testmerge250_dili3_rf) <- rownames(merged250set)
testmerge250_dili3_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili3_rf2) <- "CAM_ID"
testmerge250_dili3_rf2$Dili_3 <- as.character(testmerge250_dili3_rf[Test$CAM_ID])

testmerge250_dili3_rf2$Dili_3 <- sub(".", "", testmerge250_dili3_rf2$Dili_3)

#write.csv(trmerge250_dili3_rf2, "Dili3_TrainingPrediction_merge250rf.csv", quote=FALSE)
#


testmerge500_dili3_svm <- predict(svm500.3, merged500set)
names(testmerge500_dili3_svm) <- rownames(merged500set)
testmerge500_dili3_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili3_svm2) <- "CAM_ID"
testmerge500_dili3_svm2$Dili_3 <- as.character(testmerge500_dili3_svm[Test$CAM_ID])

testmerge500_dili3_svm2$Dili_3 <- sub(".", "", testmerge500_dili3_svm2$Dili_3)

#write.csv(testmerge500_dili3_svm2, "Dili3_TrainingPrediction_merge500.csv", quote=FALSE)
#
testmerge1000_dili3_svm <- predict(svm1000.3, merged1000set)
names(testmerge1000_dili3_svm) <- rownames(merged1000set)
testmerge1000_dili3_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili3_svm2) <- "CAM_ID"
testmerge1000_dili3_svm2$Dili_3 <- as.character(testmerge1000_dili3_svm[Test$CAM_ID])

testmerge1000_dili3_svm2$Dili_3 <- sub(".", "", testmerge1000_dili3_svm2$Dili_3)

#write.csv(trmerge1000_dili3_svm2, "Dili3_TrainingPrediction_merge1000.csv", quote=FALSE)
####
testmerge100_dili5_rf <- predict(rf100.5g, merged100set)
names(testmerge100_dili5_rf) <- rownames(merged100set)
testmerge100_dili5_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili5_rf2) <- "CAM_ID"
testmerge100_dili5_rf2$Dili_5 <- as.character(testmerge100_dili5_rf[Test$CAM_ID])

testmerge100_dili5_rf2$Dili_5 <- sub(".", "", testmerge100_dili5_rf2$Dili_5)

#write.csv(testmerge100_dili5_rf2, "Dili5_TrainingPrediction_merge100.csv", quote=FALSE)
#
testmerge250_dili5_rf <- predict(rf250.5, merged250set)
names(testmerge250_dili5_rf) <- rownames(merged250set)
testmerge250_dili5_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili5_rf2) <- "CAM_ID"
testmerge250_dili5_rf2$Dili_5 <- as.character(testmerge250_dili5_rf[Test$CAM_ID])

testmerge250_dili5_rf2$Dili_5 <- sub(".", "", testmerge250_dili5_rf2$Dili_5)

#write.csv(trmerge250_dili5_svm2, "Dili5_TrainingPrediction_merge250.csv", quote=FALSE)
#
#testmerge250_dili5_svm <- predict(svm250.5, merged250set)
#names(testmerge250_dili5_svm) <- rownames(merged250set)
#testmerge250_dili5_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
#colnames(testmerge250_dili5_svm2) <- "CAM_ID"
#testmerge250_dili5_svm2$Dili_5 <- as.character(testmerge250_dili5_svm[Test$CAM_ID])

#testmerge250_dili5_svm2$Dili_5 <- sub(".", "", testmerge250_dili5_svm2$Dili_5)

#write.csv(trmerge250_dili5_svm2, "Dili5_TrainingPrediction_merge250.csv", quote=FALSE)
#
testmerge500_dili5_svm <- predict(svm500.5, merged500set)
names(testmerge500_dili5_svm) <- rownames(merged500set)
testmerge500_dili5_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili5_svm2) <- "CAM_ID"
testmerge500_dili5_svm2$Dili_5 <- as.character(testmerge500_dili5_svm[Test$CAM_ID])

testmerge500_dili5_svm2$Dili_5 <- sub(".", "", testmerge500_dili5_svm2$Dili_5)

#write.csv(trmerge500_dili5_svm2, "Dili5_TrainingPrediction_merge500.csv", quote=FALSE)
#
testmerge1000_dili5_svm <- predict(svm1000.5, merged1000set)
names(testmerge1000_dili5_svm) <- rownames(merged1000set)
testmerge1000_dili5_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili5_svm2) <- "CAM_ID"
testmerge1000_dili5_svm2$Dili_5 <- as.character(testmerge1000_dili5_svm[Test$CAM_ID])

testmerge1000_dili5_svm2$Dili_5 <- sub(".", "", testmerge1000_dili5_svm2$Dili_5)

#write.csv(trmerge1000_dili5_svm2, "Dili5_TrainingPrediction_merge1000.csv", quote=FALSE)
####
testmerge100_dili6_rf <- predict(rf100.6g, merged100set)
names(testmerge100_dili6_rf) <- rownames(merged100set)
testmerge100_dili6_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili6_rf2) <- "CAM_ID"
testmerge100_dili6_rf2$Dili_6 <- as.character(testmerge100_dili6_rf[Test$CAM_ID])

testmerge100_dili6_rf2$Dili_6 <- sub(".", "", testmerge100_dili6_rf2$Dili_6)

#write.csv(trmerge100_dili6_svm2, "Dili6_TrainingPrediction_merge100.csv", quote=FALSE)
#
testmerge250_dili6_rf <- predict(rf250.6, merged250set)
names(testmerge250_dili6_rf) <- rownames(merged250set)
testmerge250_dili6_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili6_rf2) <- "CAM_ID"
testmerge250_dili6_rf2$Dili_6 <- as.character(testmerge250_dili6_rf[Test$CAM_ID])

testmerge250_dili6_rf2$Dili_6 <- sub(".", "", testmerge250_dili6_rf2$Dili_6)

#write.csv(trmerge250_dili6_rf2, "Dili6_TrainingPrediction_merge250.csv", quote=FALSE)
#
#testmerge250_dili6_svm <- predict(svm250.6, merged250set)
#names(testmerge250_dili6_svm) <- rownames(merged250set)
#testmerge250_dili6_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
#colnames(testmerge250_dili6_svm2) <- "CAM_ID"
#testmerge250_dili6_svm2$Dili_6 <- as.character(testmerge250_dili6_svm[Test$CAM_ID])

#testmerge250_dili6_svm2$Dili_6 <- sub(".", "", testmerge250_dili6_svm2$Dili_6)

#write.csv(trmerge250_dili6_svm2, "Dili6_TrainingPrediction_merge250.csv", quote=FALSE)
#
testmerge500_dili6_svm <- predict(svm500.6, merged500set)
names(testmerge500_dili6_svm) <- rownames(merged500set)
testmerge500_dili6_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili6_svm2) <- "CAM_ID"
testmerge500_dili6_svm2$Dili_6 <- as.character(testmerge500_dili6_svm[Test$CAM_ID])

testmerge500_dili6_svm2$Dili_6 <- sub(".", "", testmerge500_dili6_svm2$Dili_6)

#write.csv(trmerge500_dili6_svm2, "Dili6_TrainingPrediction_merge500.csv", quote=FALSE)
#
testmerge1000_dili6_svm <- predict(svm1000.6, merged1000set)
names(testmerge1000_dili6_svm) <- rownames(merged1000set)
testmerge1000_dili6_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili6_svm2) <- "CAM_ID"
testmerge1000_dili6_svm2$Dili_6 <- as.character(testmerge1000_dili6_svm[Test$CAM_ID])

testmerge1000_dili6_svm2$Dili_6 <- sub(".", "", testmerge1000_dili6_svm2$Dili_6)

#write.csv(trmerge1000_dili6_svm2, "Dili6_TrainingPrediction_merge1000.csv", quote=FALSE)

pred100 <- merge(testmerge100_dili1_svm2, testmerge100_dili3_rf2, by="CAM_ID")
pred100 <- merge(pred100, testmerge100_dili5_rf2, by="CAM_ID")
pred100 <- merge(pred100, testmerge100_dili6_rf2, by="CAM_ID")
colnames(pred100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(pred100, "../output_files/merge100_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

pred250 <- merge(testmerge250_dili1_svm2, testmerge100_dili3_rf2, by="CAM_ID")
pred250 <- merge(pred250, testmerge250_dili5_rf2, by="CAM_ID")
pred250 <- merge(pred250, testmerge250_dili6_rf2, by="CAM_ID")
colnames(pred250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(pred250, "../output_files/merge250_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

pred500 <- merge(testmerge500_dili1_svm2, testmerge500_dili3_svm2, by="CAM_ID")
pred500 <- merge(pred500, testmerge500_dili5_svm2, by="CAM_ID")
pred500 <- merge(pred500, testmerge500_dili6_svm2, by="CAM_ID")
colnames(pred500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(pred500, "../output_files/merge500_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

pred1000 <- merge(testmerge1000_dili1_nb2, testmerge1000_dili3_svm2, by="CAM_ID")
pred1000 <- merge(pred1000, testmerge500_dili5_svm2, by="CAM_ID")
pred1000 <- merge(pred1000, testmerge500_dili6_svm2, by="CAM_ID")
colnames(pred1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(pred1000, "../output_files/merge1000_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)


##MCC values

dili1_100svm_mcc <- as.data.frame(svm100.1$pred)
dili1_100svm_mcc[,1] <- as.character(dili1_100svm_mcc[,1])
dili1_100svm_mcc[,2] <- as.character(dili1_100svm_mcc[,2])
dili1_100svm_mcc[dili1_100svm_mcc=="x1"] <- 1
dili1_100svm_mcc[dili1_100svm_mcc=="x0"] <- 0

dili1_100svm<-dili1_100svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_100svm) <- c("5_fold_CV","DILI1")

dili3_100rf_mcc <- as.data.frame(rf100.3g$pred)
dili3_100rf_mcc[,1] <- as.character(dili3_100rf_mcc[,1])
dili3_100rf_mcc[,2] <- as.character(dili3_100rf_mcc[,2])
dili3_100rf_mcc[dili3_100rf_mcc=="x1"] <- 1
dili3_100rf_mcc[dili3_100rf_mcc=="x0"] <- 0

dili3_100rf<-dili3_100rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_100rf) <- c("5_fold_CV","DILI3")

#dili3_100svm_mcc <- as.data.frame(svm100.3$pred)
#dili3_100svm_mcc[,1] <- as.character(dili3_100svm_mcc[,1])
#dili3_100svm_mcc[,2] <- as.character(dili3_100svm_mcc[,2])
#dili3_100svm_mcc[dili3_100svm_mcc=="x1"] <- 1
#dili3_100svm_mcc[dili3_100svm_mcc=="x0"] <- 0
#dili3_100svm<-dili3_100svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
#  summarise(mcc=mccr(obs,pred))
#colnames(dili3_100svm) <- c("5_fold_CV","DILI3")

dili5_100rf_mcc <- as.data.frame(rf100.5g$pred)
dili5_100rf_mcc[,1] <- as.character(dili5_100rf_mcc[,1])
dili5_100rf_mcc[,2] <- as.character(dili5_100rf_mcc[,2])
dili5_100rf_mcc[dili5_100rf_mcc=="x1"] <- 1
dili5_100rf_mcc[dili5_100rf_mcc=="x0"] <- 0

dili5_100rf<-dili5_100rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_100rf) <- c("5_fold_CV","DILI5")

dili6_100rf_mcc <- as.data.frame(rf100.6g$pred)
dili6_100rf_mcc[,1] <- as.character(dili6_100rf_mcc[,1])
dili6_100rf_mcc[,2] <- as.character(dili6_100rf_mcc[,2])
dili6_100rf_mcc[dili6_100rf_mcc=="x1"] <- 1
dili6_100rf_mcc[dili6_100rf_mcc=="x0"] <- 0

dili6_100rf<-dili6_100rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_100rf) <- c("5_fold_CV","DILI6")

mcc_100merge <- merge(dili1_100svm, dili3_100rf, by="5_fold_CV")
mcc_100merge <- merge(mcc_100merge, dili5_100rf, by="5_fold_CV")
mcc_100merge <- merge(mcc_100merge, dili6_100rf, by="5_fold_CV")
write.csv(mcc_100merge, "../output_files/merge100_p3p4p5p6p7p8-crossvalidation-camda2020-UND_FEB2021.csv", quote=FALSE, row.names=FALSE)

##

dili1_250svm_mcc <- as.data.frame(svm250.1$pred)
dili1_250svm_mcc[,1] <- as.character(dili1_250svm_mcc[,1])
dili1_250svm_mcc[,2] <- as.character(dili1_250svm_mcc[,2])
dili1_250svm_mcc[dili1_250svm_mcc=="x1"] <- 1
dili1_250svm_mcc[dili1_250svm_mcc=="x0"] <- 0

dili1_250svm<-dili1_250svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_250svm) <- c("5_fold_CV","DILI1")

dili3_250rf_mcc <- as.data.frame(rf250.3$pred)
dili3_250rf_mcc[,1] <- as.character(dili3_250rf_mcc[,1])
dili3_250rf_mcc[,2] <- as.character(dili3_250rf_mcc[,2])
dili3_250rf_mcc[dili3_250rf_mcc=="x1"] <- 1
dili3_250rf_mcc[dili3_250rf_mcc=="x0"] <- 0

dili3_250rf<-dili3_250rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_250rf) <- c("5_fold_CV","DILI3")

#dili3_250svm_mcc <- as.data.frame(svm250.3$pred)
#dili3_250svm_mcc[,1] <- as.character(dili3_250svm_mcc[,1])
#dili3_250svm_mcc[,2] <- as.character(dili3_250svm_mcc[,2])
#dili3_250svm_mcc[dili3_250svm_mcc=="x1"] <- 1
#dili3_250svm_mcc[dili3_250svm_mcc=="x0"] <- 0
#dili3_250svm<-dili3_250svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
#  summarise(mcc=mccr(obs,pred))
#colnames(dili3_250svm) <- c("5_fold_CV","DILI3")

dili5_250rf_mcc <- as.data.frame(rf250.5$pred)
dili5_250rf_mcc[,1] <- as.character(dili5_250rf_mcc[,1])
dili5_250rf_mcc[,2] <- as.character(dili5_250rf_mcc[,2])
dili5_250rf_mcc[dili5_250rf_mcc=="x1"] <- 1
dili5_250rf_mcc[dili5_250rf_mcc=="x0"] <- 0

dili5_250rf<-dili5_250rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_250rf) <- c("5_fold_CV","DILI5")

dili6_250rf_mcc <- as.data.frame(rf250.6$pred)
dili6_250rf_mcc[,1] <- as.character(dili6_250rf_mcc[,1])
dili6_250rf_mcc[,2] <- as.character(dili6_250rf_mcc[,2])
dili6_250rf_mcc[dili6_250rf_mcc=="x1"] <- 1
dili6_250rf_mcc[dili6_250rf_mcc=="x0"] <- 0

dili6_250rf<-dili6_250rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_250rf) <- c("5_fold_CV","DILI6")

mcc_250merge <- merge(dili1_250svm, dili3_250rf, by="5_fold_CV")
mcc_250merge <- merge(mcc_250merge, dili5_250rf, by="5_fold_CV")
mcc_250merge <- merge(mcc_250merge, dili6_250rf, by="5_fold_CV")
write.csv(mcc_250merge, "../output_files/merge250_p3p4p5p6p7p8-crossvalidation-camda2020-UND_FEB2021.csv", quote=FALSE, row.names=FALSE)

##

dili1_500svm_mcc <- as.data.frame(svm500.1$pred)
dili1_500svm_mcc[,1] <- as.character(dili1_500svm_mcc[,1])
dili1_500svm_mcc[,2] <- as.character(dili1_500svm_mcc[,2])
dili1_500svm_mcc[dili1_500svm_mcc=="x1"] <- 1
dili1_500svm_mcc[dili1_500svm_mcc=="x0"] <- 0

dili1_500svm<-dili1_500svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_500svm) <- c("5_fold_CV","DILI1")

dili3_500svm_mcc <- as.data.frame(svm500.3$pred)
dili3_500svm_mcc[,1] <- as.character(dili3_500svm_mcc[,1])
dili3_500svm_mcc[,2] <- as.character(dili3_500svm_mcc[,2])
dili3_500svm_mcc[dili3_500svm_mcc=="x1"] <- 1
dili3_500svm_mcc[dili3_500svm_mcc=="x0"] <- 0

dili3_500svm<-dili3_500svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_500svm) <- c("5_fold_CV","DILI3")

dili5_500svm_mcc <- as.data.frame(svm500.5$pred)
dili5_500svm_mcc[,1] <- as.character(dili5_500svm_mcc[,1])
dili5_500svm_mcc[,2] <- as.character(dili5_500svm_mcc[,2])
dili5_500svm_mcc[dili5_500svm_mcc=="x1"] <- 1
dili5_500svm_mcc[dili5_500svm_mcc=="x0"] <- 0

dili5_500svm<-dili5_500svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_500svm) <- c("5_fold_CV","DILI5")

dili6_500svm_mcc <- as.data.frame(svm500.6$pred)
dili6_500svm_mcc[,1] <- as.character(dili6_500svm_mcc[,1])
dili6_500svm_mcc[,2] <- as.character(dili6_500svm_mcc[,2])
dili6_500svm_mcc[dili6_500svm_mcc=="x1"] <- 1
dili6_500svm_mcc[dili6_500svm_mcc=="x0"] <- 0

dili6_500svm<-dili6_500svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_500svm) <- c("5_fold_CV","DILI6")

mcc_500merge <- merge(dili1_500svm, dili3_500svm, by="5_fold_CV")
mcc_500merge <- merge(mcc_500merge, dili5_500svm, by="5_fold_CV")
mcc_500merge <- merge(mcc_500merge, dili6_500svm, by="5_fold_CV")
write.csv(mcc_500merge, "../output_files/merge500_p3p4p5p6p7p8-crossvalidation-camda2020-UND_FEB2021.csv", quote=FALSE, row.names=FALSE)

##

dili1_1000nb_mcc <- as.data.frame(NB1000.1$pred)
dili1_1000nb_mcc[,1] <- as.character(dili1_1000nb_mcc[,1])
dili1_1000nb_mcc[,2] <- as.character(dili1_1000nb_mcc[,2])
dili1_1000nb_mcc[dili1_1000nb_mcc=="x1"] <- 1
dili1_1000nb_mcc[dili1_1000nb_mcc=="x0"] <- 0

dili1_1000nb<-dili1_1000nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_1000nb) <- c("5_fold_CV","DILI1")

dili3_1000svm_mcc <- as.data.frame(svm1000.3$pred)
dili3_1000svm_mcc[,1] <- as.character(dili3_1000svm_mcc[,1])
dili3_1000svm_mcc[,2] <- as.character(dili3_1000svm_mcc[,2])
dili3_1000svm_mcc[dili3_1000svm_mcc=="x1"] <- 1
dili3_1000svm_mcc[dili3_1000svm_mcc=="x0"] <- 0

dili3_1000svm<-dili3_1000svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_1000svm) <- c("5_fold_CV","DILI3")

dili5_1000svm_mcc <- as.data.frame(svm1000.5$pred)
dili5_1000svm_mcc[,1] <- as.character(dili5_1000svm_mcc[,1])
dili5_1000svm_mcc[,2] <- as.character(dili5_1000svm_mcc[,2])
dili5_1000svm_mcc[dili5_1000svm_mcc=="x1"] <- 1
dili5_1000svm_mcc[dili5_1000svm_mcc=="x0"] <- 0

dili5_1000svm<-dili5_1000svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_1000svm) <- c("5_fold_CV","DILI5")

dili6_1000svm_mcc <- as.data.frame(svm1000.6$pred)
dili6_1000svm_mcc[,1] <- as.character(dili6_1000svm_mcc[,1])
dili6_1000svm_mcc[,2] <- as.character(dili6_1000svm_mcc[,2])
dili6_1000svm_mcc[dili6_1000svm_mcc=="x1"] <- 1
dili6_1000svm_mcc[dili6_1000svm_mcc=="x0"] <- 0

dili6_1000svm<-dili6_1000svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_1000svm) <- c("5_fold_CV","DILI6")

mcc_1000merge <- merge(dili1_1000nb, dili3_1000svm, by="5_fold_CV")
mcc_1000merge <- merge(mcc_1000merge, dili5_1000svm, by="5_fold_CV")
mcc_1000merge <- merge(mcc_1000merge, dili6_1000svm, by="5_fold_CV")
write.csv(mcc_1000merge, "../output_files/merge1000_p3p4p5p6p7p8-crossvalidation-camda2020-UND_FEB2021.csv", quote=FALSE, row.names=FALSE)

#####
##Predictions for all models
######
Test <- targets[which(targets$Training_Validation==""),]

modelsum <- read.csv("../output_files/ModelSummary_v2.csv")
modelsum100 <- subset(modelsum, Gene.Size=="100")
modelsum250 <- subset(modelsum, Gene.Size=="250")
modelsum500 <- subset(modelsum, Gene.Size=="500")
modelsum1000 <- subset(modelsum, Gene.Size=="1000")


#p100 <- function(x){
#  z <- x
#  predict(z, merged100set)
#  names(z) <- rownames(merged100set)
#  y <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
#  colnames(y) <- "CAM_ID"
#  w <- merge(y, z, by="CAM_ID")
#  w[,2] <- sub(".","", w[,2])
#}

#n <- p100(svm100.1)

#n <- p100(split(modelsum100, modelsum100$Model))


#lapply(split(modelsum, modelsum$Model),
#       function(x)write.table(x, file=paste(x$Model[1],".txt", sep="\t")))


testmerge100_dili1_svm <- predict(svm100.1, merged100set)
names(testmerge100_dili1_svm) <- rownames(merged100set)
testmerge100_dili1_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili1_svm2) <- "CAM_ID"
testmerge100_dili1_svm2$Dili_1 <- as.character(testmerge100_dili1_svm[Test$CAM_ID])
testmerge100_dili1_svm2$Dili_1 <- sub(".", "", testmerge100_dili1_svm2$Dili_1)

testmerge250_dili1_svm <- predict(svm250.1, merged250set)
names(testmerge250_dili1_svm) <- rownames(merged250set)
testmerge250_dili1_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili1_svm2) <- "CAM_ID"
testmerge250_dili1_svm2$Dili_1 <- as.character(testmerge250_dili1_svm[Test$CAM_ID])
testmerge250_dili1_svm2$Dili_1 <- sub(".", "", testmerge250_dili1_svm2$Dili_1)

testmerge500_dili1_svm <- predict(svm500.1, merged500set)
names(testmerge500_dili1_svm) <- rownames(merged500set)
testmerge500_dili1_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili1_svm2) <- "CAM_ID"
testmerge500_dili1_svm2$Dili_1 <- as.character(testmerge500_dili1_svm[Test$CAM_ID])
testmerge500_dili1_svm2$Dili_1 <- sub(".", "", testmerge500_dili1_svm2$Dili_1)

testmerge1000_dili1_svm <- predict(svm1000.1, merged1000set)
names(testmerge1000_dili1_svm) <- rownames(merged1000set)
testmerge1000_dili1_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili1_svm2) <- "CAM_ID"
testmerge1000_dili1_svm2$Dili_1 <- as.character(testmerge1000_dili1_svm[Test$CAM_ID])
testmerge1000_dili1_svm2$Dili_1 <- sub(".", "", testmerge1000_dili1_svm2$Dili_1)

testmerge100_dili3_svm <- predict(svm100.3, merged100set)
names(testmerge100_dili3_svm) <- rownames(merged100set)
testmerge100_dili3_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili3_svm2) <- "CAM_ID"
testmerge100_dili3_svm2$Dili_3 <- as.character(testmerge100_dili3_svm[Test$CAM_ID])
testmerge100_dili3_svm2$Dili_3 <- sub(".", "", testmerge100_dili3_svm2$Dili_3)

testmerge250_dili3_svm <- predict(svm250.3, merged250set)
names(testmerge250_dili3_svm) <- rownames(merged250set)
testmerge250_dili3_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili3_svm2) <- "CAM_ID"
testmerge250_dili3_svm2$Dili_3 <- as.character(testmerge250_dili3_svm[Test$CAM_ID])
testmerge250_dili3_svm2$Dili_3 <- sub(".", "", testmerge250_dili3_svm2$Dili_3)

testmerge500_dili3_svm <- predict(svm500.3, merged500set)
names(testmerge500_dili3_svm) <- rownames(merged500set)
testmerge500_dili3_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili3_svm2) <- "CAM_ID"
testmerge500_dili3_svm2$Dili_3 <- as.character(testmerge500_dili3_svm[Test$CAM_ID])
testmerge500_dili3_svm2$Dili_3 <- sub(".", "", testmerge500_dili3_svm2$Dili_3)

testmerge1000_dili3_svm <- predict(svm1000.3, merged1000set)
names(testmerge1000_dili3_svm) <- rownames(merged1000set)
testmerge1000_dili3_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili3_svm2) <- "CAM_ID"
testmerge1000_dili3_svm2$Dili_3 <- as.character(testmerge1000_dili3_svm[Test$CAM_ID])
testmerge1000_dili3_svm2$Dili_3 <- sub(".", "", testmerge1000_dili3_svm2$Dili_3)

testmerge100_dili5_svm <- predict(svm100.5, merged100set)
names(testmerge100_dili5_svm) <- rownames(merged100set)
testmerge100_dili5_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili5_svm2) <- "CAM_ID"
testmerge100_dili5_svm2$Dili_5 <- as.character(testmerge100_dili5_svm[Test$CAM_ID])
testmerge100_dili5_svm2$Dili_5 <- sub(".", "", testmerge100_dili5_svm2$Dili_5)

testmerge250_dili5_svm <- predict(svm250.5, merged250set)
names(testmerge250_dili5_svm) <- rownames(merged250set)
testmerge250_dili5_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili5_svm2) <- "CAM_ID"
testmerge250_dili5_svm2$Dili_5 <- as.character(testmerge250_dili5_svm[Test$CAM_ID])
testmerge250_dili5_svm2$Dili_5 <- sub(".", "", testmerge250_dili5_svm2$Dili_5)

testmerge500_dili5_svm <- predict(svm500.5, merged500set)
names(testmerge500_dili5_svm) <- rownames(merged500set)
testmerge500_dili5_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili5_svm2) <- "CAM_ID"
testmerge500_dili5_svm2$Dili_5 <- as.character(testmerge500_dili5_svm[Test$CAM_ID])
testmerge500_dili5_svm2$Dili_5 <- sub(".", "", testmerge500_dili5_svm2$Dili_5)

testmerge1000_dili5_svm <- predict(svm1000.5, merged1000set)
names(testmerge1000_dili5_svm) <- rownames(merged1000set)
testmerge1000_dili5_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili5_svm2) <- "CAM_ID"
testmerge1000_dili5_svm2$Dili_5 <- as.character(testmerge1000_dili5_svm[Test$CAM_ID])
testmerge1000_dili5_svm2$Dili_5 <- sub(".", "", testmerge1000_dili5_svm2$Dili_5)

testmerge100_dili6_svm <- predict(svm100.6, merged100set)
names(testmerge100_dili6_svm) <- rownames(merged100set)
testmerge100_dili6_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili6_svm2) <- "CAM_ID"
testmerge100_dili6_svm2$Dili_6 <- as.character(testmerge100_dili6_svm[Test$CAM_ID])
testmerge100_dili6_svm2$Dili_6 <- sub(".", "", testmerge100_dili6_svm2$Dili_6)

testmerge250_dili6_svm <- predict(svm250.6, merged250set)
names(testmerge250_dili6_svm) <- rownames(merged250set)
testmerge250_dili6_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili6_svm2) <- "CAM_ID"
testmerge250_dili6_svm2$Dili_6 <- as.character(testmerge250_dili6_svm[Test$CAM_ID])
testmerge250_dili6_svm2$Dili_6 <- sub(".", "", testmerge250_dili6_svm2$Dili_6)

testmerge500_dili6_svm <- predict(svm500.6, merged500set)
names(testmerge500_dili6_svm) <- rownames(merged500set)
testmerge500_dili6_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili6_svm2) <- "CAM_ID"
testmerge500_dili6_svm2$Dili_6 <- as.character(testmerge500_dili6_svm[Test$CAM_ID])
testmerge500_dili6_svm2$Dili_6 <- sub(".", "", testmerge500_dili6_svm2$Dili_6)

testmerge1000_dili6_svm <- predict(svm1000.6, merged1000set)
names(testmerge1000_dili6_svm) <- rownames(merged1000set)
testmerge1000_dili6_svm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili6_svm2) <- "CAM_ID"
testmerge1000_dili6_svm2$Dili_6 <- as.character(testmerge1000_dili6_svm[Test$CAM_ID])
testmerge1000_dili6_svm2$Dili_6 <- sub(".", "", testmerge1000_dili6_svm2$Dili_6)

testmerge100_dili1_rf <- predict(rf100.1g, merged100set)
names(testmerge100_dili1_rf) <- rownames(merged100set)
testmerge100_dili1_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili1_rf2) <- "CAM_ID"
testmerge100_dili1_rf2$Dili_1 <- as.character(testmerge100_dili1_rf[Test$CAM_ID])
testmerge100_dili1_rf2$Dili_1 <- sub(".", "", testmerge100_dili1_rf2$Dili_1)

testmerge250_dili1_rf <- predict(rf250.1, merged250set)
names(testmerge250_dili1_rf) <- rownames(merged250set)
testmerge250_dili1_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili1_rf2) <- "CAM_ID"
testmerge250_dili1_rf2$Dili_1 <- as.character(testmerge250_dili1_rf[Test$CAM_ID])
testmerge250_dili1_rf2$Dili_1 <- sub(".", "", testmerge250_dili1_rf2$Dili_1)

testmerge500_dili1_rf <- predict(rf500.1, merged500set)
names(testmerge500_dili1_rf) <- rownames(merged500set)
testmerge500_dili1_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili1_rf2) <- "CAM_ID"
testmerge500_dili1_rf2$Dili_1 <- as.character(testmerge500_dili1_rf[Test$CAM_ID])
testmerge500_dili1_rf2$Dili_1 <- sub(".", "", testmerge500_dili1_rf2$Dili_1)

testmerge1000_dili1_rf <- predict(rf1000.1, merged1000set)
names(testmerge1000_dili1_rf) <- rownames(merged1000set)
testmerge1000_dili1_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili1_rf2) <- "CAM_ID"
testmerge1000_dili1_rf2$Dili_1 <- as.character(testmerge1000_dili1_rf[Test$CAM_ID])
testmerge1000_dili1_rf2$Dili_1 <- sub(".", "", testmerge1000_dili1_rf2$Dili_1)

testmerge100_dili3_rf <- predict(rf100.3g, merged100set)
names(testmerge100_dili3_rf) <- rownames(merged100set)
testmerge100_dili3_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili3_rf2) <- "CAM_ID"
testmerge100_dili3_rf2$Dili_3 <- as.character(testmerge100_dili3_rf[Test$CAM_ID])
testmerge100_dili3_rf2$Dili_3 <- sub(".", "", testmerge100_dili3_rf2$Dili_3)

testmerge250_dili3_rf <- predict(rf250.3, merged250set)
names(testmerge250_dili3_rf) <- rownames(merged250set)
testmerge250_dili3_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili3_rf2) <- "CAM_ID"
testmerge250_dili3_rf2$Dili_3 <- as.character(testmerge250_dili3_rf[Test$CAM_ID])
testmerge250_dili3_rf2$Dili_3 <- sub(".", "", testmerge250_dili3_rf2$Dili_3)

testmerge500_dili3_rf <- predict(rf500.3, merged500set)
names(testmerge500_dili3_rf) <- rownames(merged500set)
testmerge500_dili3_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili3_rf2) <- "CAM_ID"
testmerge500_dili3_rf2$Dili_3 <- as.character(testmerge500_dili3_rf[Test$CAM_ID])
testmerge500_dili3_rf2$Dili_3 <- sub(".", "", testmerge500_dili3_rf2$Dili_3)

testmerge1000_dili3_rf <- predict(rf1000.3, merged1000set)
names(testmerge1000_dili3_rf) <- rownames(merged1000set)
testmerge1000_dili3_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili3_rf2) <- "CAM_ID"
testmerge1000_dili3_rf2$Dili_3 <- as.character(testmerge1000_dili3_rf[Test$CAM_ID])
testmerge1000_dili3_rf2$Dili_3 <- sub(".", "", testmerge1000_dili3_rf2$Dili_3)

testmerge100_dili5_rf <- predict(rf100.5g, merged100set)
names(testmerge100_dili5_rf) <- rownames(merged100set)
testmerge100_dili5_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili5_rf2) <- "CAM_ID"
testmerge100_dili5_rf2$Dili_5 <- as.character(testmerge100_dili5_rf[Test$CAM_ID])
testmerge100_dili5_rf2$Dili_5 <- sub(".", "", testmerge100_dili5_rf2$Dili_5)

testmerge250_dili5_rf <- predict(rf250.5, merged250set)
names(testmerge250_dili5_rf) <- rownames(merged250set)
testmerge250_dili5_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili5_rf2) <- "CAM_ID"
testmerge250_dili5_rf2$Dili_5 <- as.character(testmerge250_dili5_rf[Test$CAM_ID])
testmerge250_dili5_rf2$Dili_5 <- sub(".", "", testmerge250_dili5_rf2$Dili_5)

testmerge500_dili5_rf <- predict(rf500.5, merged500set)
names(testmerge500_dili5_rf) <- rownames(merged500set)
testmerge500_dili5_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili5_rf2) <- "CAM_ID"
testmerge500_dili5_rf2$Dili_5 <- as.character(testmerge500_dili5_rf[Test$CAM_ID])
testmerge500_dili5_rf2$Dili_5 <- sub(".", "", testmerge500_dili5_rf2$Dili_5)

testmerge1000_dili5_rf <- predict(rf1000.5, merged1000set)
names(testmerge1000_dili5_rf) <- rownames(merged1000set)
testmerge1000_dili5_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili5_rf2) <- "CAM_ID"
testmerge1000_dili5_rf2$Dili_5 <- as.character(testmerge1000_dili5_rf[Test$CAM_ID])
testmerge1000_dili5_rf2$Dili_5 <- sub(".", "", testmerge1000_dili5_rf2$Dili_5)

testmerge100_dili6_rf <- predict(rf100.6g, merged100set)
names(testmerge100_dili6_rf) <- rownames(merged100set)
testmerge100_dili6_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili6_rf2) <- "CAM_ID"
testmerge100_dili6_rf2$Dili_6 <- as.character(testmerge100_dili6_rf[Test$CAM_ID])
testmerge100_dili6_rf2$Dili_6 <- sub(".", "", testmerge100_dili6_rf2$Dili_6)

testmerge250_dili6_rf <- predict(rf250.6, merged250set)
names(testmerge250_dili6_rf) <- rownames(merged250set)
testmerge250_dili6_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili6_rf2) <- "CAM_ID"
testmerge250_dili6_rf2$Dili_6 <- as.character(testmerge250_dili6_rf[Test$CAM_ID])
testmerge250_dili6_rf2$Dili_6 <- sub(".", "", testmerge250_dili6_rf2$Dili_6)

testmerge500_dili6_rf <- predict(rf500.6, merged500set)
names(testmerge500_dili6_rf) <- rownames(merged500set)
testmerge500_dili6_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili6_rf2) <- "CAM_ID"
testmerge500_dili6_rf2$Dili_6 <- as.character(testmerge500_dili6_rf[Test$CAM_ID])
testmerge500_dili6_rf2$Dili_6 <- sub(".", "", testmerge500_dili6_rf2$Dili_6)

testmerge1000_dili6_rf <- predict(rf1000.6, merged1000set)
names(testmerge1000_dili6_rf) <- rownames(merged1000set)
testmerge1000_dili6_rf2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili6_rf2) <- "CAM_ID"
testmerge1000_dili6_rf2$Dili_6 <- as.character(testmerge1000_dili6_rf[Test$CAM_ID])
testmerge1000_dili6_rf2$Dili_6 <- sub(".", "", testmerge1000_dili6_rf2$Dili_6)

testmerge100_dili1_glm <- predict(glm100.1, merged100set)
names(testmerge100_dili1_glm) <- rownames(merged100set)
testmerge100_dili1_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili1_glm2) <- "CAM_ID"
testmerge100_dili1_glm2$Dili_1 <- as.character(testmerge100_dili1_glm[Test$CAM_ID])
testmerge100_dili1_glm2$Dili_1 <- sub(".", "", testmerge100_dili1_glm2$Dili_1)

testmerge250_dili1_glm <- predict(glm250.1, merged250set)
names(testmerge250_dili1_glm) <- rownames(merged250set)
testmerge250_dili1_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili1_glm2) <- "CAM_ID"
testmerge250_dili1_glm2$Dili_1 <- as.character(testmerge250_dili1_glm[Test$CAM_ID])
testmerge250_dili1_glm2$Dili_1 <- sub(".", "", testmerge250_dili1_glm2$Dili_1)

testmerge500_dili1_glm <- predict(glm500.1, merged500set)
names(testmerge500_dili1_glm) <- rownames(merged500set)
testmerge500_dili1_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili1_glm2) <- "CAM_ID"
testmerge500_dili1_glm2$Dili_1 <- as.character(testmerge500_dili1_glm[Test$CAM_ID])
testmerge500_dili1_glm2$Dili_1 <- sub(".", "", testmerge500_dili1_glm2$Dili_1)

testmerge1000_dili1_glm <- predict(glm1000.1, merged1000set)
names(testmerge1000_dili1_glm) <- rownames(merged1000set)
testmerge1000_dili1_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili1_glm2) <- "CAM_ID"
testmerge1000_dili1_glm2$Dili_1 <- as.character(testmerge1000_dili1_glm[Test$CAM_ID])
testmerge1000_dili1_glm2$Dili_1 <- sub(".", "", testmerge1000_dili1_glm2$Dili_1)

testmerge100_dili3_glm <- predict(glm100.3, merged100set)
names(testmerge100_dili3_glm) <- rownames(merged100set)
testmerge100_dili3_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili3_glm2) <- "CAM_ID"
testmerge100_dili3_glm2$Dili_3 <- as.character(testmerge100_dili3_glm[Test$CAM_ID])
testmerge100_dili3_glm2$Dili_3 <- sub(".", "", testmerge100_dili3_glm2$Dili_3)

testmerge250_dili3_glm <- predict(glm250.3, merged250set)
names(testmerge250_dili3_glm) <- rownames(merged250set)
testmerge250_dili3_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili3_glm2) <- "CAM_ID"
testmerge250_dili3_glm2$Dili_3 <- as.character(testmerge250_dili3_glm[Test$CAM_ID])
testmerge250_dili3_glm2$Dili_3 <- sub(".", "", testmerge250_dili3_glm2$Dili_3)

testmerge500_dili3_glm <- predict(glm500.3, merged500set)
names(testmerge500_dili3_glm) <- rownames(merged500set)
testmerge500_dili3_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili3_glm2) <- "CAM_ID"
testmerge500_dili3_glm2$Dili_3 <- as.character(testmerge500_dili3_glm[Test$CAM_ID])
testmerge500_dili3_glm2$Dili_3 <- sub(".", "", testmerge500_dili3_glm2$Dili_3)

testmerge1000_dili3_glm <- predict(glm1000.3, merged1000set)
names(testmerge1000_dili3_glm) <- rownames(merged1000set)
testmerge1000_dili3_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili3_glm2) <- "CAM_ID"
testmerge1000_dili3_glm2$Dili_3 <- as.character(testmerge1000_dili3_glm[Test$CAM_ID])
testmerge1000_dili3_glm2$Dili_3 <- sub(".", "", testmerge1000_dili3_glm2$Dili_3)

testmerge100_dili5_glm <- predict(glm100.5, merged100set)
names(testmerge100_dili5_glm) <- rownames(merged100set)
testmerge100_dili5_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili5_glm2) <- "CAM_ID"
testmerge100_dili5_glm2$Dili_5 <- as.character(testmerge100_dili5_glm[Test$CAM_ID])
testmerge100_dili5_glm2$Dili_5 <- sub(".", "", testmerge100_dili5_glm2$Dili_5)

testmerge250_dili5_glm <- predict(glm250.5, merged250set)
names(testmerge250_dili5_glm) <- rownames(merged250set)
testmerge250_dili5_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili5_glm2) <- "CAM_ID"
testmerge250_dili5_glm2$Dili_5 <- as.character(testmerge250_dili5_glm[Test$CAM_ID])
testmerge250_dili5_glm2$Dili_5 <- sub(".", "", testmerge250_dili5_glm2$Dili_5)

testmerge500_dili5_glm <- predict(glm500.5, merged500set)
names(testmerge500_dili5_glm) <- rownames(merged500set)
testmerge500_dili5_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili5_glm2) <- "CAM_ID"
testmerge500_dili5_glm2$Dili_5 <- as.character(testmerge500_dili5_glm[Test$CAM_ID])
testmerge500_dili5_glm2$Dili_5 <- sub(".", "", testmerge500_dili5_glm2$Dili_5)

testmerge1000_dili5_glm <- predict(glm1000.5, merged1000set)
names(testmerge1000_dili5_glm) <- rownames(merged1000set)
testmerge1000_dili5_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili5_glm2) <- "CAM_ID"
testmerge1000_dili5_glm2$Dili_5 <- as.character(testmerge1000_dili5_glm[Test$CAM_ID])
testmerge1000_dili5_glm2$Dili_5 <- sub(".", "", testmerge1000_dili5_glm2$Dili_5)

testmerge100_dili6_glm <- predict(glm100.6, merged100set)
names(testmerge100_dili6_glm) <- rownames(merged100set)
testmerge100_dili6_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili6_glm2) <- "CAM_ID"
testmerge100_dili6_glm2$Dili_6 <- as.character(testmerge100_dili6_glm[Test$CAM_ID])
testmerge100_dili6_glm2$Dili_6 <- sub(".", "", testmerge100_dili6_glm2$Dili_6)

testmerge250_dili6_glm <- predict(glm250.6, merged250set)
names(testmerge250_dili6_glm) <- rownames(merged250set)
testmerge250_dili6_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili6_glm2) <- "CAM_ID"
testmerge250_dili6_glm2$Dili_6 <- as.character(testmerge250_dili6_glm[Test$CAM_ID])
testmerge250_dili6_glm2$Dili_6 <- sub(".", "", testmerge250_dili6_glm2$Dili_6)

testmerge500_dili6_glm <- predict(glm500.6, merged500set)
names(testmerge500_dili6_glm) <- rownames(merged500set)
testmerge500_dili6_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili6_glm2) <- "CAM_ID"
testmerge500_dili6_glm2$Dili_6 <- as.character(testmerge500_dili6_glm[Test$CAM_ID])
testmerge500_dili6_glm2$Dili_6 <- sub(".", "", testmerge500_dili6_glm2$Dili_6)

testmerge1000_dili6_glm <- predict(glm1000.6, merged1000set)
names(testmerge1000_dili6_glm) <- rownames(merged1000set)
testmerge1000_dili6_glm2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili6_glm2) <- "CAM_ID"
testmerge1000_dili6_glm2$Dili_6 <- as.character(testmerge1000_dili6_glm[Test$CAM_ID])
testmerge1000_dili6_glm2$Dili_6 <- sub(".", "", testmerge1000_dili6_glm2$Dili_6)

testmerge100_dili1_rpart <- predict(rpart100.1, merged100set)
names(testmerge100_dili1_rpart) <- rownames(merged100set)
testmerge100_dili1_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili1_rpart2) <- "CAM_ID"
testmerge100_dili1_rpart2$Dili_1 <- as.character(testmerge100_dili1_rpart[Test$CAM_ID])
testmerge100_dili1_rpart2$Dili_1 <- sub(".", "", testmerge100_dili1_rpart2$Dili_1)

testmerge250_dili1_rpart <- predict(rpart250.1, merged250set)
names(testmerge250_dili1_rpart) <- rownames(merged250set)
testmerge250_dili1_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili1_rpart2) <- "CAM_ID"
testmerge250_dili1_rpart2$Dili_1 <- as.character(testmerge250_dili1_rpart[Test$CAM_ID])
testmerge250_dili1_rpart2$Dili_1 <- sub(".", "", testmerge250_dili1_rpart2$Dili_1)

testmerge500_dili1_rpart <- predict(rpart500.1, merged500set)
names(testmerge500_dili1_rpart) <- rownames(merged500set)
testmerge500_dili1_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili1_rpart2) <- "CAM_ID"
testmerge500_dili1_rpart2$Dili_1 <- as.character(testmerge500_dili1_rpart[Test$CAM_ID])
testmerge500_dili1_rpart2$Dili_1 <- sub(".", "", testmerge500_dili1_rpart2$Dili_1)

testmerge1000_dili1_rpart <- predict(rpart1000.1, merged1000set)
names(testmerge1000_dili1_rpart) <- rownames(merged1000set)
testmerge1000_dili1_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili1_rpart2) <- "CAM_ID"
testmerge1000_dili1_rpart2$Dili_1 <- as.character(testmerge1000_dili1_rpart[Test$CAM_ID])
testmerge1000_dili1_rpart2$Dili_1 <- sub(".", "", testmerge1000_dili1_rpart2$Dili_1)

testmerge100_dili3_rpart <- predict(rpart100.3, merged100set)
names(testmerge100_dili3_rpart) <- rownames(merged100set)
testmerge100_dili3_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili3_rpart2) <- "CAM_ID"
testmerge100_dili3_rpart2$Dili_3 <- as.character(testmerge100_dili3_rpart[Test$CAM_ID])
testmerge100_dili3_rpart2$Dili_3 <- sub(".", "", testmerge100_dili3_rpart2$Dili_3)

testmerge250_dili3_rpart <- predict(rpart250.3, merged250set)
names(testmerge250_dili3_rpart) <- rownames(merged250set)
testmerge250_dili3_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili3_rpart2) <- "CAM_ID"
testmerge250_dili3_rpart2$Dili_3 <- as.character(testmerge250_dili3_rpart[Test$CAM_ID])
testmerge250_dili3_rpart2$Dili_3 <- sub(".", "", testmerge250_dili3_rpart2$Dili_3)

testmerge500_dili3_rpart <- predict(rpart500.3, merged500set)
names(testmerge500_dili3_rpart) <- rownames(merged500set)
testmerge500_dili3_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili3_rpart2) <- "CAM_ID"
testmerge500_dili3_rpart2$Dili_3 <- as.character(testmerge500_dili3_rpart[Test$CAM_ID])
testmerge500_dili3_rpart2$Dili_3 <- sub(".", "", testmerge500_dili3_rpart2$Dili_3)

testmerge1000_dili3_rpart <- predict(rpart1000.3, merged1000set)
names(testmerge1000_dili3_rpart) <- rownames(merged1000set)
testmerge1000_dili3_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili3_rpart2) <- "CAM_ID"
testmerge1000_dili3_rpart2$Dili_3 <- as.character(testmerge1000_dili3_rpart[Test$CAM_ID])
testmerge1000_dili3_rpart2$Dili_3 <- sub(".", "", testmerge1000_dili3_rpart2$Dili_3)

testmerge100_dili5_rpart <- predict(rpart100.5, merged100set)
names(testmerge100_dili5_rpart) <- rownames(merged100set)
testmerge100_dili5_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili5_rpart2) <- "CAM_ID"
testmerge100_dili5_rpart2$Dili_5 <- as.character(testmerge100_dili5_rpart[Test$CAM_ID])
testmerge100_dili5_rpart2$Dili_5 <- sub(".", "", testmerge100_dili5_rpart2$Dili_5)

testmerge250_dili5_rpart <- predict(rpart250.5, merged250set)
names(testmerge250_dili5_rpart) <- rownames(merged250set)
testmerge250_dili5_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili5_rpart2) <- "CAM_ID"
testmerge250_dili5_rpart2$Dili_5 <- as.character(testmerge250_dili5_rpart[Test$CAM_ID])
testmerge250_dili5_rpart2$Dili_5 <- sub(".", "", testmerge250_dili5_rpart2$Dili_5)

testmerge500_dili5_rpart <- predict(rpart500.5, merged500set)
names(testmerge500_dili5_rpart) <- rownames(merged500set)
testmerge500_dili5_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili5_rpart2) <- "CAM_ID"
testmerge500_dili5_rpart2$Dili_5 <- as.character(testmerge500_dili5_rpart[Test$CAM_ID])
testmerge500_dili5_rpart2$Dili_5 <- sub(".", "", testmerge500_dili5_rpart2$Dili_5)

testmerge1000_dili5_rpart <- predict(rpart1000.5, merged1000set)
names(testmerge1000_dili5_rpart) <- rownames(merged1000set)
testmerge1000_dili5_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili5_rpart2) <- "CAM_ID"
testmerge1000_dili5_rpart2$Dili_5 <- as.character(testmerge1000_dili5_rpart[Test$CAM_ID])
testmerge1000_dili5_rpart2$Dili_5 <- sub(".", "", testmerge1000_dili5_rpart2$Dili_5)

testmerge100_dili6_rpart <- predict(rpart100.6, merged100set)
names(testmerge100_dili6_rpart) <- rownames(merged100set)
testmerge100_dili6_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili6_rpart2) <- "CAM_ID"
testmerge100_dili6_rpart2$Dili_6 <- as.character(testmerge100_dili6_rpart[Test$CAM_ID])
testmerge100_dili6_rpart2$Dili_6 <- sub(".", "", testmerge100_dili6_rpart2$Dili_6)

testmerge250_dili6_rpart <- predict(rpart250.6, merged250set)
names(testmerge250_dili6_rpart) <- rownames(merged250set)
testmerge250_dili6_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili6_rpart2) <- "CAM_ID"
testmerge250_dili6_rpart2$Dili_6 <- as.character(testmerge250_dili6_rpart[Test$CAM_ID])
testmerge250_dili6_rpart2$Dili_6 <- sub(".", "", testmerge250_dili6_rpart2$Dili_6)

testmerge500_dili6_rpart <- predict(rpart500.6, merged500set)
names(testmerge500_dili6_rpart) <- rownames(merged500set)
testmerge500_dili6_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili6_rpart2) <- "CAM_ID"
testmerge500_dili6_rpart2$Dili_6 <- as.character(testmerge500_dili6_rpart[Test$CAM_ID])
testmerge500_dili6_rpart2$Dili_6 <- sub(".", "", testmerge500_dili6_rpart2$Dili_6)

testmerge1000_dili6_rpart <- predict(rpart1000.6, merged1000set)
names(testmerge1000_dili6_rpart) <- rownames(merged1000set)
testmerge1000_dili6_rpart2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili6_rpart2) <- "CAM_ID"
testmerge1000_dili6_rpart2$Dili_6 <- as.character(testmerge1000_dili6_rpart[Test$CAM_ID])
testmerge1000_dili6_rpart2$Dili_6 <- sub(".", "", testmerge1000_dili6_rpart2$Dili_6)

testmerge100_dili1_nb <- predict(NB100.1, merged100set)
names(testmerge100_dili1_nb) <- rownames(merged100set)
testmerge100_dili1_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili1_nb2) <- "CAM_ID"
testmerge100_dili1_nb2$Dili_1 <- as.character(testmerge100_dili1_nb[Test$CAM_ID])
testmerge100_dili1_nb2$Dili_1 <- sub(".", "", testmerge100_dili1_nb2$Dili_1)

testmerge250_dili1_nb <- predict(NB250.1, merged250set)
names(testmerge250_dili1_nb) <- rownames(merged250set)
testmerge250_dili1_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili1_nb2) <- "CAM_ID"
testmerge250_dili1_nb2$Dili_1 <- as.character(testmerge250_dili1_nb[Test$CAM_ID])
testmerge250_dili1_nb2$Dili_1 <- sub(".", "", testmerge250_dili1_nb2$Dili_1)

testmerge500_dili1_nb <- predict(NB500.1, merged500set)
names(testmerge500_dili1_nb) <- rownames(merged500set)
testmerge500_dili1_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili1_nb2) <- "CAM_ID"
testmerge500_dili1_nb2$Dili_1 <- as.character(testmerge500_dili1_nb[Test$CAM_ID])
testmerge500_dili1_nb2$Dili_1 <- sub(".", "", testmerge500_dili1_nb2$Dili_1)

testmerge1000_dili1_nb <- predict(NB1000.1, merged1000set)
names(testmerge1000_dili1_nb) <- rownames(merged1000set)
testmerge1000_dili1_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili1_nb2) <- "CAM_ID"
testmerge1000_dili1_nb2$Dili_1 <- as.character(testmerge1000_dili1_nb[Test$CAM_ID])
testmerge1000_dili1_nb2$Dili_1 <- sub(".", "", testmerge1000_dili1_nb2$Dili_1)

testmerge100_dili3_nb <- predict(NB100.3, merged100set)
names(testmerge100_dili3_nb) <- rownames(merged100set)
testmerge100_dili3_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili3_nb2) <- "CAM_ID"
testmerge100_dili3_nb2$Dili_3 <- as.character(testmerge100_dili3_nb[Test$CAM_ID])
testmerge100_dili3_nb2$Dili_3 <- sub(".", "", testmerge100_dili3_nb2$Dili_3)

testmerge250_dili3_nb <- predict(NB250.3, merged250set)
names(testmerge250_dili3_nb) <- rownames(merged250set)
testmerge250_dili3_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili3_nb2) <- "CAM_ID"
testmerge250_dili3_nb2$Dili_3 <- as.character(testmerge250_dili3_nb[Test$CAM_ID])
testmerge250_dili3_nb2$Dili_3 <- sub(".", "", testmerge250_dili3_nb2$Dili_3)

testmerge500_dili3_nb <- predict(NB500.3, merged500set)
names(testmerge500_dili3_nb) <- rownames(merged500set)
testmerge500_dili3_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili3_nb2) <- "CAM_ID"
testmerge500_dili3_nb2$Dili_3 <- as.character(testmerge500_dili3_nb[Test$CAM_ID])
testmerge500_dili3_nb2$Dili_3 <- sub(".", "", testmerge500_dili3_nb2$Dili_3)

testmerge1000_dili3_nb <- predict(NB1000.3, merged1000set)
names(testmerge1000_dili3_nb) <- rownames(merged1000set)
testmerge1000_dili3_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili3_nb2) <- "CAM_ID"
testmerge1000_dili3_nb2$Dili_3 <- as.character(testmerge1000_dili3_nb[Test$CAM_ID])
testmerge1000_dili3_nb2$Dili_3 <- sub(".", "", testmerge1000_dili3_nb2$Dili_3)

testmerge100_dili5_nb <- predict(NB100.5, merged100set)
names(testmerge100_dili5_nb) <- rownames(merged100set)
testmerge100_dili5_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili5_nb2) <- "CAM_ID"
testmerge100_dili5_nb2$Dili_5 <- as.character(testmerge100_dili5_nb[Test$CAM_ID])
testmerge100_dili5_nb2$Dili_5 <- sub(".", "", testmerge100_dili5_nb2$Dili_5)

testmerge250_dili5_nb <- predict(NB250.5, merged250set)
names(testmerge250_dili5_nb) <- rownames(merged250set)
testmerge250_dili5_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili5_nb2) <- "CAM_ID"
testmerge250_dili5_nb2$Dili_5 <- as.character(testmerge250_dili5_nb[Test$CAM_ID])
testmerge250_dili5_nb2$Dili_5 <- sub(".", "", testmerge250_dili5_nb2$Dili_5)

testmerge500_dili5_nb <- predict(NB500.5, merged500set)
names(testmerge500_dili5_nb) <- rownames(merged500set)
testmerge500_dili5_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili5_nb2) <- "CAM_ID"
testmerge500_dili5_nb2$Dili_5 <- as.character(testmerge500_dili5_nb[Test$CAM_ID])
testmerge500_dili5_nb2$Dili_5 <- sub(".", "", testmerge500_dili5_nb2$Dili_5)

testmerge1000_dili5_nb <- predict(NB1000.5, merged1000set)
names(testmerge1000_dili5_nb) <- rownames(merged1000set)
testmerge1000_dili5_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili5_nb2) <- "CAM_ID"
testmerge1000_dili5_nb2$Dili_5 <- as.character(testmerge1000_dili5_nb[Test$CAM_ID])
testmerge1000_dili5_nb2$Dili_5 <- sub(".", "", testmerge1000_dili5_nb2$Dili_5)

testmerge100_dili6_nb <- predict(NB100.6, merged100set)
names(testmerge100_dili6_nb) <- rownames(merged100set)
testmerge100_dili6_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge100_dili6_nb2) <- "CAM_ID"
testmerge100_dili6_nb2$Dili_6 <- as.character(testmerge100_dili6_nb[Test$CAM_ID])
testmerge100_dili6_nb2$Dili_6 <- sub(".", "", testmerge100_dili6_nb2$Dili_6)

testmerge250_dili6_nb <- predict(NB250.6, merged250set)
names(testmerge250_dili6_nb) <- rownames(merged250set)
testmerge250_dili6_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge250_dili6_nb2) <- "CAM_ID"
testmerge250_dili6_nb2$Dili_6 <- as.character(testmerge250_dili6_nb[Test$CAM_ID])
testmerge250_dili6_nb2$Dili_6 <- sub(".", "", testmerge250_dili6_nb2$Dili_6)

testmerge500_dili6_nb <- predict(NB500.6, merged500set)
names(testmerge500_dili6_nb) <- rownames(merged500set)
testmerge500_dili6_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge500_dili6_nb2) <- "CAM_ID"
testmerge500_dili6_nb2$Dili_6 <- as.character(testmerge500_dili6_nb[Test$CAM_ID])
testmerge500_dili6_nb2$Dili_6 <- sub(".", "", testmerge500_dili6_nb2$Dili_6)

testmerge1000_dili6_nb <- predict(NB1000.6, merged1000set)
names(testmerge1000_dili6_nb) <- rownames(merged1000set)
testmerge1000_dili6_nb2 <- as.data.frame(Test$CAM_ID, stringsAsFactors=FALSE)
colnames(testmerge1000_dili6_nb2) <- "CAM_ID"
testmerge1000_dili6_nb2$Dili_6 <- as.character(testmerge1000_dili6_nb[Test$CAM_ID])
testmerge1000_dili6_nb2$Dili_6 <- sub(".", "", testmerge1000_dili6_nb2$Dili_6)

predsvm100 <- merge(testmerge100_dili1_svm2, testmerge100_dili3_svm2, by="CAM_ID")
predsvm100 <- merge(predsvm100, testmerge100_dili5_svm2, by="CAM_ID")
predsvm100 <- merge(predsvm100, testmerge100_dili6_svm2, by="CAM_ID")
colnames(predsvm100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predsvm100, "../output_files/merge100svm_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predsvm250 <- merge(testmerge250_dili1_svm2, testmerge250_dili3_svm2, by="CAM_ID")
predsvm250 <- merge(predsvm250, testmerge250_dili5_svm2, by="CAM_ID")
predsvm250 <- merge(predsvm250, testmerge250_dili6_svm2, by="CAM_ID")
colnames(predsvm250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predsvm250, "../output_files/merge250svm_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predsvm500 <- merge(testmerge500_dili1_svm2, testmerge500_dili3_svm2, by="CAM_ID")
predsvm500 <- merge(predsvm500, testmerge500_dili5_svm2, by="CAM_ID")
predsvm500 <- merge(predsvm500, testmerge500_dili6_svm2, by="CAM_ID")
colnames(predsvm500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predsvm500, "../output_files/merge500svm_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predsvm1000 <- merge(testmerge1000_dili1_svm2, testmerge1000_dili3_svm2, by="CAM_ID")
predsvm1000 <- merge(predsvm1000, testmerge1000_dili5_svm2, by="CAM_ID")
predsvm1000 <- merge(predsvm1000, testmerge1000_dili6_svm2, by="CAM_ID")
colnames(predsvm1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predsvm1000, "../output_files/merge1000svm_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predrf100 <- merge(testmerge100_dili1_rf2, testmerge100_dili3_rf2, by="CAM_ID")
predrf100 <- merge(predrf100, testmerge100_dili5_rf2, by="CAM_ID")
predrf100 <- merge(predrf100, testmerge100_dili6_rf2, by="CAM_ID")
colnames(predrf100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predrf100, "../output_files/merge100rf_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predrf250 <- merge(testmerge250_dili1_rf2, testmerge250_dili3_rf2, by="CAM_ID")
predrf250 <- merge(predrf250, testmerge250_dili5_rf2, by="CAM_ID")
predrf250 <- merge(predrf250, testmerge250_dili6_rf2, by="CAM_ID")
colnames(predrf250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predrf250, "../output_files/merge250rf_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predrf500 <- merge(testmerge500_dili1_rf2, testmerge500_dili3_rf2, by="CAM_ID")
predrf500 <- merge(predrf500, testmerge500_dili5_rf2, by="CAM_ID")
predrf500 <- merge(predrf500, testmerge500_dili6_rf2, by="CAM_ID")
colnames(predrf500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predrf500, "../output_files/merge500rf_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predrf1000 <- merge(testmerge1000_dili1_rf2, testmerge1000_dili3_rf2, by="CAM_ID")
predrf1000 <- merge(predrf1000, testmerge1000_dili5_rf2, by="CAM_ID")
predrf1000 <- merge(predrf1000, testmerge1000_dili6_rf2, by="CAM_ID")
colnames(predrf1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predrf1000, "../output_files/merge1000rf_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predglm100 <- merge(testmerge100_dili1_glm2, testmerge100_dili3_glm2, by="CAM_ID")
predglm100 <- merge(predglm100, testmerge100_dili5_glm2, by="CAM_ID")
predglm100 <- merge(predglm100, testmerge100_dili6_glm2, by="CAM_ID")
colnames(predglm100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predglm100, "../output_files/merge100glm_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predglm250 <- merge(testmerge250_dili1_glm2, testmerge250_dili3_glm2, by="CAM_ID")
predglm250 <- merge(predglm250, testmerge250_dili5_glm2, by="CAM_ID")
predglm250 <- merge(predglm250, testmerge250_dili6_glm2, by="CAM_ID")
colnames(predglm250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predglm250, "../output_files/merge250glm_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predglm500 <- merge(testmerge500_dili1_glm2, testmerge500_dili3_glm2, by="CAM_ID")
predglm500 <- merge(predglm500, testmerge500_dili5_glm2, by="CAM_ID")
predglm500 <- merge(predglm500, testmerge500_dili6_glm2, by="CAM_ID")
colnames(predglm500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predglm500, "../output_files/merge500glm_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predglm1000 <- merge(testmerge1000_dili1_glm2, testmerge1000_dili3_glm2, by="CAM_ID")
predglm1000 <- merge(predglm1000, testmerge1000_dili5_glm2, by="CAM_ID")
predglm1000 <- merge(predglm1000, testmerge1000_dili6_glm2, by="CAM_ID")
colnames(predglm1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predglm1000, "../output_files/merge1000glm_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predrpart100 <- merge(testmerge100_dili1_rpart2, testmerge100_dili3_rpart2, by="CAM_ID")
predrpart100 <- merge(predrpart100, testmerge100_dili5_rpart2, by="CAM_ID")
predrpart100 <- merge(predrpart100, testmerge100_dili6_rpart2, by="CAM_ID")
colnames(predrpart100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predrpart100, "../output_files/merge100rpart_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predrpart250 <- merge(testmerge250_dili1_rpart2, testmerge250_dili3_rpart2, by="CAM_ID")
predrpart250 <- merge(predrpart250, testmerge250_dili5_rpart2, by="CAM_ID")
predrpart250 <- merge(predrpart250, testmerge250_dili6_rpart2, by="CAM_ID")
colnames(predrpart250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predrpart250, "../output_files/merge250rpart_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predrpart500 <- merge(testmerge500_dili1_rpart2, testmerge500_dili3_rpart2, by="CAM_ID")
predrpart500 <- merge(predrpart500, testmerge500_dili5_rpart2, by="CAM_ID")
predrpart500 <- merge(predrpart500, testmerge500_dili6_rpart2, by="CAM_ID")
colnames(predrpart500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predrpart500, "../output_files/merge500rpart_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

predrpart1000 <- merge(testmerge1000_dili1_rpart2, testmerge1000_dili3_rpart2, by="CAM_ID")
predrpart1000 <- merge(predrpart1000, testmerge1000_dili5_rpart2, by="CAM_ID")
predrpart1000 <- merge(predrpart1000, testmerge1000_dili6_rpart2, by="CAM_ID")
colnames(predrpart1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(predrpart1000, "../output_files/merge1000rpart_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

prednb100 <- merge(testmerge100_dili1_nb2, testmerge100_dili3_nb2, by="CAM_ID")
prednb100 <- merge(prednb100, testmerge100_dili5_nb2, by="CAM_ID")
prednb100 <- merge(prednb100, testmerge100_dili6_nb2, by="CAM_ID")
colnames(prednb100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(prednb100, "../output_files/merge100nb_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

prednb250 <- merge(testmerge250_dili1_nb2, testmerge250_dili3_nb2, by="CAM_ID")
prednb250 <- merge(prednb250, testmerge250_dili5_nb2, by="CAM_ID")
prednb250 <- merge(prednb250, testmerge250_dili6_nb2, by="CAM_ID")
colnames(prednb250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(prednb250, "../output_files/merge250nb_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

prednb500 <- merge(testmerge500_dili1_nb2, testmerge500_dili3_nb2, by="CAM_ID")
prednb500 <- merge(prednb500, testmerge500_dili5_nb2, by="CAM_ID")
prednb500 <- merge(prednb500, testmerge500_dili6_nb2, by="CAM_ID")
colnames(prednb500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(prednb500, "../output_files/merge500nb_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

prednb1000 <- merge(testmerge1000_dili1_nb2, testmerge1000_dili3_nb2, by="CAM_ID")
prednb1000 <- merge(prednb1000, testmerge1000_dili5_nb2, by="CAM_ID")
prednb1000 <- merge(prednb1000, testmerge1000_dili6_nb2, by="CAM_ID")
colnames(prednb1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(prednb1000, "../output_files/merge1000nb_p3p4p5p6p7p8-predictions-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

#######
##MCC All Models
######

dili1_100svm_mcc <- as.data.frame(svm100.1$pred)
dili1_100svm_mcc[,1] <- as.character(dili1_100svm_mcc[,1])
dili1_100svm_mcc[,2] <- as.character(dili1_100svm_mcc[,2])
dili1_100svm_mcc[dili1_100svm_mcc=="x1"] <- 1
dili1_100svm_mcc[dili1_100svm_mcc=="x0"] <- 0

dili1_100svm<- dili1_100svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_100svm) <- c("5_fold_CV","DILI1")

dili1_250svm_mcc <- as.data.frame(svm250.1$pred)
dili1_250svm_mcc[,1] <- as.character(dili1_250svm_mcc[,1])
dili1_250svm_mcc[,2] <- as.character(dili1_250svm_mcc[,2])
dili1_250svm_mcc[dili1_250svm_mcc=="x1"] <- 1
dili1_250svm_mcc[dili1_250svm_mcc=="x0"] <- 0

dili1_250svm<-dili1_250svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_250svm) <- c("5_fold_CV","DILI1")

dili1_500svm_mcc <- as.data.frame(svm500.1$pred)
dili1_500svm_mcc[,1] <- as.character(dili1_500svm_mcc[,1])
dili1_500svm_mcc[,2] <- as.character(dili1_500svm_mcc[,2])
dili1_500svm_mcc[dili1_500svm_mcc=="x1"] <- 1
dili1_500svm_mcc[dili1_500svm_mcc=="x0"] <- 0

dili1_500svm<-dili1_500svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_500svm) <- c("5_fold_CV","DILI1")

dili1_1000svm_mcc <- as.data.frame(svm1000.1$pred)
dili1_1000svm_mcc[,1] <- as.character(dili1_1000svm_mcc[,1])
dili1_1000svm_mcc[,2] <- as.character(dili1_1000svm_mcc[,2])
dili1_1000svm_mcc[dili1_1000svm_mcc=="x1"] <- 1
dili1_1000svm_mcc[dili1_1000svm_mcc=="x0"] <- 0

dili1_1000svm<-dili1_1000svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_1000svm) <- c("5_fold_CV","DILI1")

dili3_100svm_mcc <- as.data.frame(svm100.3$pred)
dili3_100svm_mcc[,1] <- as.character(dili3_100svm_mcc[,1])
dili3_100svm_mcc[,2] <- as.character(dili3_100svm_mcc[,2])
dili3_100svm_mcc[dili3_100svm_mcc=="x1"] <- 1
dili3_100svm_mcc[dili3_100svm_mcc=="x0"] <- 0

dili3_100svm<-dili3_100svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_100svm) <- c("5_fold_CV","DILI3")

dili3_250svm_mcc <- as.data.frame(svm250.3$pred)
dili3_250svm_mcc[,1] <- as.character(dili3_250svm_mcc[,1])
dili3_250svm_mcc[,2] <- as.character(dili3_250svm_mcc[,2])
dili3_250svm_mcc[dili3_250svm_mcc=="x1"] <- 1
dili3_250svm_mcc[dili3_250svm_mcc=="x0"] <- 0

dili3_250svm<-dili3_250svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_250svm) <- c("5_fold_CV","DILI3")

dili3_500svm_mcc <- as.data.frame(svm500.3$pred)
dili3_500svm_mcc[,1] <- as.character(dili3_500svm_mcc[,1])
dili3_500svm_mcc[,2] <- as.character(dili3_500svm_mcc[,2])
dili3_500svm_mcc[dili3_500svm_mcc=="x1"] <- 1
dili3_500svm_mcc[dili3_500svm_mcc=="x0"] <- 0

dili3_500svm<-dili3_500svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_500svm) <- c("5_fold_CV","DILI3")

dili3_1000svm_mcc <- as.data.frame(svm1000.3$pred)
dili3_1000svm_mcc[,1] <- as.character(dili3_1000svm_mcc[,1])
dili3_1000svm_mcc[,2] <- as.character(dili3_1000svm_mcc[,2])
dili3_1000svm_mcc[dili3_1000svm_mcc=="x1"] <- 1
dili3_1000svm_mcc[dili3_1000svm_mcc=="x0"] <- 0

dili3_1000svm<-dili3_1000svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_1000svm) <- c("5_fold_CV","DILI3")

dili5_100svm_mcc <- as.data.frame(svm100.5$pred)
dili5_100svm_mcc[,1] <- as.character(dili5_100svm_mcc[,1])
dili5_100svm_mcc[,2] <- as.character(dili5_100svm_mcc[,2])
dili5_100svm_mcc[dili5_100svm_mcc=="x1"] <- 1
dili5_100svm_mcc[dili5_100svm_mcc=="x0"] <- 0

dili5_100svm<-dili5_100svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_100svm) <- c("5_fold_CV","DILI5")

dili5_250svm_mcc <- as.data.frame(svm250.5$pred)
dili5_250svm_mcc[,1] <- as.character(dili5_250svm_mcc[,1])
dili5_250svm_mcc[,2] <- as.character(dili5_250svm_mcc[,2])
dili5_250svm_mcc[dili5_250svm_mcc=="x1"] <- 1
dili5_250svm_mcc[dili5_250svm_mcc=="x0"] <- 0

dili5_250svm<-dili5_250svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_250svm) <- c("5_fold_CV","DILI5")

dili5_500svm_mcc <- as.data.frame(svm500.5$pred)
dili5_500svm_mcc[,1] <- as.character(dili5_500svm_mcc[,1])
dili5_500svm_mcc[,2] <- as.character(dili5_500svm_mcc[,2])
dili5_500svm_mcc[dili5_500svm_mcc=="x1"] <- 1
dili5_500svm_mcc[dili5_500svm_mcc=="x0"] <- 0

dili5_500svm<-dili5_500svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_500svm) <- c("5_fold_CV","DILI5")

dili5_1000svm_mcc <- as.data.frame(svm1000.5$pred)
dili5_1000svm_mcc[,1] <- as.character(dili5_1000svm_mcc[,1])
dili5_1000svm_mcc[,2] <- as.character(dili5_1000svm_mcc[,2])
dili5_1000svm_mcc[dili5_1000svm_mcc=="x1"] <- 1
dili5_1000svm_mcc[dili5_1000svm_mcc=="x0"] <- 0

dili5_1000svm<-dili5_1000svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_1000svm) <- c("5_fold_CV","DILI5")

dili6_100svm_mcc <- as.data.frame(svm100.6$pred)
dili6_100svm_mcc[,1] <- as.character(dili6_100svm_mcc[,1])
dili6_100svm_mcc[,2] <- as.character(dili6_100svm_mcc[,2])
dili6_100svm_mcc[dili6_100svm_mcc=="x1"] <- 1
dili6_100svm_mcc[dili6_100svm_mcc=="x0"] <- 0

dili6_100svm<-dili6_100svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_100svm) <- c("5_fold_CV","DILI6")

dili6_250svm_mcc <- as.data.frame(svm250.6$pred)
dili6_250svm_mcc[,1] <- as.character(dili6_250svm_mcc[,1])
dili6_250svm_mcc[,2] <- as.character(dili6_250svm_mcc[,2])
dili6_250svm_mcc[dili6_250svm_mcc=="x1"] <- 1
dili6_250svm_mcc[dili6_250svm_mcc=="x0"] <- 0

dili6_250svm<-dili6_250svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_250svm) <- c("5_fold_CV","DILI6")

dili6_500svm_mcc <- as.data.frame(svm500.6$pred)
dili6_500svm_mcc[,1] <- as.character(dili6_500svm_mcc[,1])
dili6_500svm_mcc[,2] <- as.character(dili6_500svm_mcc[,2])
dili6_500svm_mcc[dili6_500svm_mcc=="x1"] <- 1
dili6_500svm_mcc[dili6_500svm_mcc=="x0"] <- 0

dili6_500svm<-dili6_500svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_500svm) <- c("5_fold_CV","DILI6")

dili6_1000svm_mcc <- as.data.frame(svm1000.6$pred)
dili6_1000svm_mcc[,1] <- as.character(dili6_1000svm_mcc[,1])
dili6_1000svm_mcc[,2] <- as.character(dili6_1000svm_mcc[,2])
dili6_1000svm_mcc[dili6_1000svm_mcc=="x1"] <- 1
dili6_1000svm_mcc[dili6_1000svm_mcc=="x0"] <- 0

dili6_1000svm<-dili6_1000svm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_1000svm) <- c("5_fold_CV","DILI6")

dili1_100rf_mcc <- as.data.frame(rf100.1g$pred)
dili1_100rf_mcc[,1] <- as.character(dili1_100rf_mcc[,1])
dili1_100rf_mcc[,2] <- as.character(dili1_100rf_mcc[,2])
dili1_100rf_mcc[dili1_100rf_mcc=="x1"] <- 1
dili1_100rf_mcc[dili1_100rf_mcc=="x0"] <- 0

dili1_100rf<-dili1_100rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_100rf) <- c("5_fold_CV","DILI1")

dili1_250rf_mcc <- as.data.frame(rf250.1$pred)
dili1_250rf_mcc[,1] <- as.character(dili1_250rf_mcc[,1])
dili1_250rf_mcc[,2] <- as.character(dili1_250rf_mcc[,2])
dili1_250rf_mcc[dili1_250rf_mcc=="x1"] <- 1
dili1_250rf_mcc[dili1_250rf_mcc=="x0"] <- 0

dili1_250rf<-dili1_250rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_250rf) <- c("5_fold_CV","DILI1")

dili1_500rf_mcc <- as.data.frame(rf500.1$pred)
dili1_500rf_mcc[,1] <- as.character(dili1_500rf_mcc[,1])
dili1_500rf_mcc[,2] <- as.character(dili1_500rf_mcc[,2])
dili1_500rf_mcc[dili1_500rf_mcc=="x1"] <- 1
dili1_500rf_mcc[dili1_500rf_mcc=="x0"] <- 0

dili1_500rf<-dili1_500rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_500rf) <- c("5_fold_CV","DILI1")

dili1_1000rf_mcc <- as.data.frame(rf1000.1$pred)
dili1_1000rf_mcc[,1] <- as.character(dili1_1000rf_mcc[,1])
dili1_1000rf_mcc[,2] <- as.character(dili1_1000rf_mcc[,2])
dili1_1000rf_mcc[dili1_1000rf_mcc=="x1"] <- 1
dili1_1000rf_mcc[dili1_1000rf_mcc=="x0"] <- 0

dili1_1000rf<-dili1_1000rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_1000rf) <- c("5_fold_CV","DILI1")

dili3_100rf_mcc <- as.data.frame(rf100.3g$pred)
dili3_100rf_mcc[,1] <- as.character(dili3_100rf_mcc[,1])
dili3_100rf_mcc[,2] <- as.character(dili3_100rf_mcc[,2])
dili3_100rf_mcc[dili3_100rf_mcc=="x1"] <- 1
dili3_100rf_mcc[dili3_100rf_mcc=="x0"] <- 0

dili3_100rf<-dili3_100rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_100rf) <- c("5_fold_CV","DILI3")

dili3_250rf_mcc <- as.data.frame(rf250.3$pred)
dili3_250rf_mcc[,1] <- as.character(dili3_250rf_mcc[,1])
dili3_250rf_mcc[,2] <- as.character(dili3_250rf_mcc[,2])
dili3_250rf_mcc[dili3_250rf_mcc=="x1"] <- 1
dili3_250rf_mcc[dili3_250rf_mcc=="x0"] <- 0

dili3_250rf<-dili3_250rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_250rf) <- c("5_fold_CV","DILI3")

dili3_500rf_mcc <- as.data.frame(rf500.3$pred)
dili3_500rf_mcc[,1] <- as.character(dili3_500rf_mcc[,1])
dili3_500rf_mcc[,2] <- as.character(dili3_500rf_mcc[,2])
dili3_500rf_mcc[dili3_500rf_mcc=="x1"] <- 1
dili3_500rf_mcc[dili3_500rf_mcc=="x0"] <- 0

dili3_500rf<-dili3_500rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_500rf) <- c("5_fold_CV","DILI3")

dili3_1000rf_mcc <- as.data.frame(rf1000.3$pred)
dili3_1000rf_mcc[,1] <- as.character(dili3_1000rf_mcc[,1])
dili3_1000rf_mcc[,2] <- as.character(dili3_1000rf_mcc[,2])
dili3_1000rf_mcc[dili3_1000rf_mcc=="x1"] <- 1
dili3_1000rf_mcc[dili3_1000rf_mcc=="x0"] <- 0

dili3_1000rf<-dili3_1000rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_1000rf) <- c("5_fold_CV","DILI3")

dili5_100rf_mcc <- as.data.frame(rf100.5g$pred)
dili5_100rf_mcc[,1] <- as.character(dili5_100rf_mcc[,1])
dili5_100rf_mcc[,2] <- as.character(dili5_100rf_mcc[,2])
dili5_100rf_mcc[dili5_100rf_mcc=="x1"] <- 1
dili5_100rf_mcc[dili5_100rf_mcc=="x0"] <- 0

dili5_100rf<-dili5_100rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_100rf) <- c("5_fold_CV","DILI5")

dili5_250rf_mcc <- as.data.frame(rf250.5$pred)
dili5_250rf_mcc[,1] <- as.character(dili5_250rf_mcc[,1])
dili5_250rf_mcc[,2] <- as.character(dili5_250rf_mcc[,2])
dili5_250rf_mcc[dili5_250rf_mcc=="x1"] <- 1
dili5_250rf_mcc[dili5_250rf_mcc=="x0"] <- 0

dili5_250rf<-dili5_250rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_250rf) <- c("5_fold_CV","DILI5")

dili5_500rf_mcc <- as.data.frame(rf500.5$pred)
dili5_500rf_mcc[,1] <- as.character(dili5_500rf_mcc[,1])
dili5_500rf_mcc[,2] <- as.character(dili5_500rf_mcc[,2])
dili5_500rf_mcc[dili5_500rf_mcc=="x1"] <- 1
dili5_500rf_mcc[dili5_500rf_mcc=="x0"] <- 0

dili5_500rf<-dili5_500rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_500rf) <- c("5_fold_CV","DILI5")

dili5_1000rf_mcc <- as.data.frame(rf1000.5$pred)
dili5_1000rf_mcc[,1] <- as.character(dili5_1000rf_mcc[,1])
dili5_1000rf_mcc[,2] <- as.character(dili5_1000rf_mcc[,2])
dili5_1000rf_mcc[dili5_1000rf_mcc=="x1"] <- 1
dili5_1000rf_mcc[dili5_1000rf_mcc=="x0"] <- 0

dili5_1000rf<-dili5_1000rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_1000rf) <- c("5_fold_CV","DILI5")

dili6_100rf_mcc <- as.data.frame(rf100.6g$pred)
dili6_100rf_mcc[,1] <- as.character(dili6_100rf_mcc[,1])
dili6_100rf_mcc[,2] <- as.character(dili6_100rf_mcc[,2])
dili6_100rf_mcc[dili6_100rf_mcc=="x1"] <- 1
dili6_100rf_mcc[dili6_100rf_mcc=="x0"] <- 0

dili6_100rf<-dili6_100rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_100rf) <- c("5_fold_CV","DILI6")

dili6_250rf_mcc <- as.data.frame(rf250.6$pred)
dili6_250rf_mcc[,1] <- as.character(dili6_250rf_mcc[,1])
dili6_250rf_mcc[,2] <- as.character(dili6_250rf_mcc[,2])
dili6_250rf_mcc[dili6_250rf_mcc=="x1"] <- 1
dili6_250rf_mcc[dili6_250rf_mcc=="x0"] <- 0

dili6_250rf<-dili6_250rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_250rf) <- c("5_fold_CV","DILI6")

dili6_500rf_mcc <- as.data.frame(rf500.6$pred)
dili6_500rf_mcc[,1] <- as.character(dili6_500rf_mcc[,1])
dili6_500rf_mcc[,2] <- as.character(dili6_500rf_mcc[,2])
dili6_500rf_mcc[dili6_500rf_mcc=="x1"] <- 1
dili6_500rf_mcc[dili6_500rf_mcc=="x0"] <- 0

dili6_500rf<-dili6_500rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_500rf) <- c("5_fold_CV","DILI6")

dili6_1000rf_mcc <- as.data.frame(rf1000.6$pred)
dili6_1000rf_mcc[,1] <- as.character(dili6_1000rf_mcc[,1])
dili6_1000rf_mcc[,2] <- as.character(dili6_1000rf_mcc[,2])
dili6_1000rf_mcc[dili6_1000rf_mcc=="x1"] <- 1
dili6_1000rf_mcc[dili6_1000rf_mcc=="x0"] <- 0

dili6_1000rf<-dili6_1000rf_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_1000rf) <- c("5_fold_CV","DILI6")

dili1_100glm_mcc <- as.data.frame(glm100.1$pred)
dili1_100glm_mcc[,1] <- as.character(dili1_100glm_mcc[,1])
dili1_100glm_mcc[,2] <- as.character(dili1_100glm_mcc[,2])
dili1_100glm_mcc[dili1_100glm_mcc=="x1"] <- 1
dili1_100glm_mcc[dili1_100glm_mcc=="x0"] <- 0

dili1_100glm<-dili1_100glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_100glm) <- c("5_fold_CV","DILI1")

dili1_250glm_mcc <- as.data.frame(glm250.1$pred)
dili1_250glm_mcc[,1] <- as.character(dili1_250glm_mcc[,1])
dili1_250glm_mcc[,2] <- as.character(dili1_250glm_mcc[,2])
dili1_250glm_mcc[dili1_250glm_mcc=="x1"] <- 1
dili1_250glm_mcc[dili1_250glm_mcc=="x0"] <- 0

dili1_250glm<-dili1_250glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_250glm) <- c("5_fold_CV","DILI1")

dili1_500glm_mcc <- as.data.frame(glm500.1$pred)
dili1_500glm_mcc[,1] <- as.character(dili1_500glm_mcc[,1])
dili1_500glm_mcc[,2] <- as.character(dili1_500glm_mcc[,2])
dili1_500glm_mcc[dili1_500glm_mcc=="x1"] <- 1
dili1_500glm_mcc[dili1_500glm_mcc=="x0"] <- 0

dili1_500glm<-dili1_500glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_500glm) <- c("5_fold_CV","DILI1")

dili1_1000glm_mcc <- as.data.frame(glm1000.1$pred)
dili1_1000glm_mcc[,1] <- as.character(dili1_1000glm_mcc[,1])
dili1_1000glm_mcc[,2] <- as.character(dili1_1000glm_mcc[,2])
dili1_1000glm_mcc[dili1_1000glm_mcc=="x1"] <- 1
dili1_1000glm_mcc[dili1_1000glm_mcc=="x0"] <- 0

dili1_1000glm<-dili1_1000glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_1000glm) <- c("5_fold_CV","DILI1")

dili3_100glm_mcc <- as.data.frame(glm100.3$pred)
dili3_100glm_mcc[,1] <- as.character(dili3_100glm_mcc[,1])
dili3_100glm_mcc[,2] <- as.character(dili3_100glm_mcc[,2])
dili3_100glm_mcc[dili3_100glm_mcc=="x1"] <- 1
dili3_100glm_mcc[dili3_100glm_mcc=="x0"] <- 0

dili3_100glm<-dili3_100glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_100glm) <- c("5_fold_CV","DILI3")

dili3_250glm_mcc <- as.data.frame(glm250.3$pred)
dili3_250glm_mcc[,1] <- as.character(dili3_250glm_mcc[,1])
dili3_250glm_mcc[,2] <- as.character(dili3_250glm_mcc[,2])
dili3_250glm_mcc[dili3_250glm_mcc=="x1"] <- 1
dili3_250glm_mcc[dili3_250glm_mcc=="x0"] <- 0

dili3_250glm<-dili3_250glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_250glm) <- c("5_fold_CV","DILI3")

dili3_500glm_mcc <- as.data.frame(glm500.3$pred)
dili3_500glm_mcc[,1] <- as.character(dili3_500glm_mcc[,1])
dili3_500glm_mcc[,2] <- as.character(dili3_500glm_mcc[,2])
dili3_500glm_mcc[dili3_500glm_mcc=="x1"] <- 1
dili3_500glm_mcc[dili3_500glm_mcc=="x0"] <- 0

dili3_500glm<-dili3_500glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_500glm) <- c("5_fold_CV","DILI3")

dili3_1000glm_mcc <- as.data.frame(glm1000.3$pred)
dili3_1000glm_mcc[,1] <- as.character(dili3_1000glm_mcc[,1])
dili3_1000glm_mcc[,2] <- as.character(dili3_1000glm_mcc[,2])
dili3_1000glm_mcc[dili3_1000glm_mcc=="x1"] <- 1
dili3_1000glm_mcc[dili3_1000glm_mcc=="x0"] <- 0

dili3_1000glm<-dili3_1000glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_1000glm) <- c("5_fold_CV","DILI3")

dili5_100glm_mcc <- as.data.frame(glm100.5$pred)
dili5_100glm_mcc[,1] <- as.character(dili5_100glm_mcc[,1])
dili5_100glm_mcc[,2] <- as.character(dili5_100glm_mcc[,2])
dili5_100glm_mcc[dili5_100glm_mcc=="x1"] <- 1
dili5_100glm_mcc[dili5_100glm_mcc=="x0"] <- 0

dili5_100glm<-dili5_100glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_100glm) <- c("5_fold_CV","DILI5")

dili5_250glm_mcc <- as.data.frame(glm250.5$pred)
dili5_250glm_mcc[,1] <- as.character(dili5_250glm_mcc[,1])
dili5_250glm_mcc[,2] <- as.character(dili5_250glm_mcc[,2])
dili5_250glm_mcc[dili5_250glm_mcc=="x1"] <- 1
dili5_250glm_mcc[dili5_250glm_mcc=="x0"] <- 0

dili5_250glm<-dili5_250glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_250glm) <- c("5_fold_CV","DILI5")

dili5_500glm_mcc <- as.data.frame(glm500.5$pred)
dili5_500glm_mcc[,1] <- as.character(dili5_500glm_mcc[,1])
dili5_500glm_mcc[,2] <- as.character(dili5_500glm_mcc[,2])
dili5_500glm_mcc[dili5_500glm_mcc=="x1"] <- 1
dili5_500glm_mcc[dili5_500glm_mcc=="x0"] <- 0

dili5_500glm<-dili5_500glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_500glm) <- c("5_fold_CV","DILI5")

dili5_1000glm_mcc <- as.data.frame(glm1000.5$pred)
dili5_1000glm_mcc[,1] <- as.character(dili5_1000glm_mcc[,1])
dili5_1000glm_mcc[,2] <- as.character(dili5_1000glm_mcc[,2])
dili5_1000glm_mcc[dili5_1000glm_mcc=="x1"] <- 1
dili5_1000glm_mcc[dili5_1000glm_mcc=="x0"] <- 0

dili5_1000glm<-dili5_1000glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_1000glm) <- c("5_fold_CV","DILI5")

dili6_100glm_mcc <- as.data.frame(glm100.6$pred)
dili6_100glm_mcc[,1] <- as.character(dili6_100glm_mcc[,1])
dili6_100glm_mcc[,2] <- as.character(dili6_100glm_mcc[,2])
dili6_100glm_mcc[dili6_100glm_mcc=="x1"] <- 1
dili6_100glm_mcc[dili6_100glm_mcc=="x0"] <- 0

dili6_100glm<-dili6_100glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_100glm) <- c("5_fold_CV","DILI6")

dili6_250glm_mcc <- as.data.frame(glm250.6$pred)
dili6_250glm_mcc[,1] <- as.character(dili6_250glm_mcc[,1])
dili6_250glm_mcc[,2] <- as.character(dili6_250glm_mcc[,2])
dili6_250glm_mcc[dili6_250glm_mcc=="x1"] <- 1
dili6_250glm_mcc[dili6_250glm_mcc=="x0"] <- 0

dili6_250glm<-dili6_250glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_250glm) <- c("5_fold_CV","DILI6")

dili6_500glm_mcc <- as.data.frame(glm500.6$pred)
dili6_500glm_mcc[,1] <- as.character(dili6_500glm_mcc[,1])
dili6_500glm_mcc[,2] <- as.character(dili6_500glm_mcc[,2])
dili6_500glm_mcc[dili6_500glm_mcc=="x1"] <- 1
dili6_500glm_mcc[dili6_500glm_mcc=="x0"] <- 0

dili6_500glm<-dili6_500glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_500glm) <- c("5_fold_CV","DILI6")

dili6_1000glm_mcc <- as.data.frame(glm1000.6$pred)
dili6_1000glm_mcc[,1] <- as.character(dili6_1000glm_mcc[,1])
dili6_1000glm_mcc[,2] <- as.character(dili6_1000glm_mcc[,2])
dili6_1000glm_mcc[dili6_1000glm_mcc=="x1"] <- 1
dili6_1000glm_mcc[dili6_1000glm_mcc=="x0"] <- 0

dili6_1000glm<-dili6_1000glm_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_1000glm) <- c("5_fold_CV","DILI6")

dili1_100rpart_mcc <- as.data.frame(rpart100.1$pred)
dili1_100rpart_mcc[,1] <- as.character(dili1_100rpart_mcc[,1])
dili1_100rpart_mcc[,2] <- as.character(dili1_100rpart_mcc[,2])
dili1_100rpart_mcc[dili1_100rpart_mcc=="x1"] <- 1
dili1_100rpart_mcc[dili1_100rpart_mcc=="x0"] <- 0

dili1_100rpart<-dili1_100rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_100rpart) <- c("5_fold_CV","DILI1")

dili1_250rpart_mcc <- as.data.frame(rpart250.1$pred)
dili1_250rpart_mcc[,1] <- as.character(dili1_250rpart_mcc[,1])
dili1_250rpart_mcc[,2] <- as.character(dili1_250rpart_mcc[,2])
dili1_250rpart_mcc[dili1_250rpart_mcc=="x1"] <- 1
dili1_250rpart_mcc[dili1_250rpart_mcc=="x0"] <- 0

dili1_250rpart<-dili1_250rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_250rpart) <- c("5_fold_CV","DILI1")

dili1_500rpart_mcc <- as.data.frame(rpart500.1$pred)
dili1_500rpart_mcc[,1] <- as.character(dili1_500rpart_mcc[,1])
dili1_500rpart_mcc[,2] <- as.character(dili1_500rpart_mcc[,2])
dili1_500rpart_mcc[dili1_500rpart_mcc=="x1"] <- 1
dili1_500rpart_mcc[dili1_500rpart_mcc=="x0"] <- 0

dili1_500rpart<-dili1_500rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_500rpart) <- c("5_fold_CV","DILI1")

dili1_1000rpart_mcc <- as.data.frame(rpart1000.1$pred)
dili1_1000rpart_mcc[,1] <- as.character(dili1_1000rpart_mcc[,1])
dili1_1000rpart_mcc[,2] <- as.character(dili1_1000rpart_mcc[,2])
dili1_1000rpart_mcc[dili1_1000rpart_mcc=="x1"] <- 1
dili1_1000rpart_mcc[dili1_1000rpart_mcc=="x0"] <- 0

dili1_1000rpart<-dili1_1000rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_1000rpart) <- c("5_fold_CV","DILI1")

dili3_100rpart_mcc <- as.data.frame(rpart100.3$pred)
dili3_100rpart_mcc[,1] <- as.character(dili3_100rpart_mcc[,1])
dili3_100rpart_mcc[,2] <- as.character(dili3_100rpart_mcc[,2])
dili3_100rpart_mcc[dili3_100rpart_mcc=="x1"] <- 1
dili3_100rpart_mcc[dili3_100rpart_mcc=="x0"] <- 0

dili3_100rpart<-dili3_100rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_100rpart) <- c("5_fold_CV","DILI3")

dili3_250rpart_mcc <- as.data.frame(rpart250.3$pred)
dili3_250rpart_mcc[,1] <- as.character(dili3_250rpart_mcc[,1])
dili3_250rpart_mcc[,2] <- as.character(dili3_250rpart_mcc[,2])
dili3_250rpart_mcc[dili3_250rpart_mcc=="x1"] <- 1
dili3_250rpart_mcc[dili3_250rpart_mcc=="x0"] <- 0

dili3_250rpart<-dili3_250rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_250rpart) <- c("5_fold_CV","DILI3")

dili3_500rpart_mcc <- as.data.frame(rpart500.3$pred)
dili3_500rpart_mcc[,1] <- as.character(dili3_500rpart_mcc[,1])
dili3_500rpart_mcc[,2] <- as.character(dili3_500rpart_mcc[,2])
dili3_500rpart_mcc[dili3_500rpart_mcc=="x1"] <- 1
dili3_500rpart_mcc[dili3_500rpart_mcc=="x0"] <- 0

dili3_500rpart<-dili3_500rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_500rpart) <- c("5_fold_CV","DILI3")

dili3_1000rpart_mcc <- as.data.frame(rpart1000.3$pred)
dili3_1000rpart_mcc[,1] <- as.character(dili3_1000rpart_mcc[,1])
dili3_1000rpart_mcc[,2] <- as.character(dili3_1000rpart_mcc[,2])
dili3_1000rpart_mcc[dili3_1000rpart_mcc=="x1"] <- 1
dili3_1000rpart_mcc[dili3_1000rpart_mcc=="x0"] <- 0

dili3_1000rpart<-dili3_1000rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_1000rpart) <- c("5_fold_CV","DILI3")

dili5_100rpart_mcc <- as.data.frame(rpart100.5$pred)
dili5_100rpart_mcc[,1] <- as.character(dili5_100rpart_mcc[,1])
dili5_100rpart_mcc[,2] <- as.character(dili5_100rpart_mcc[,2])
dili5_100rpart_mcc[dili5_100rpart_mcc=="x1"] <- 1
dili5_100rpart_mcc[dili5_100rpart_mcc=="x0"] <- 0

dili5_100rpart<-dili5_100rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_100rpart) <- c("5_fold_CV","DILI5")

dili5_250rpart_mcc <- as.data.frame(rpart250.5$pred)
dili5_250rpart_mcc[,1] <- as.character(dili5_250rpart_mcc[,1])
dili5_250rpart_mcc[,2] <- as.character(dili5_250rpart_mcc[,2])
dili5_250rpart_mcc[dili5_250rpart_mcc=="x1"] <- 1
dili5_250rpart_mcc[dili5_250rpart_mcc=="x0"] <- 0

dili5_250rpart<-dili5_250rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_250rpart) <- c("5_fold_CV","DILI5")

dili5_500rpart_mcc <- as.data.frame(rpart500.5$pred)
dili5_500rpart_mcc[,1] <- as.character(dili5_500rpart_mcc[,1])
dili5_500rpart_mcc[,2] <- as.character(dili5_500rpart_mcc[,2])
dili5_500rpart_mcc[dili5_500rpart_mcc=="x1"] <- 1
dili5_500rpart_mcc[dili5_500rpart_mcc=="x0"] <- 0

dili5_500rpart<-dili5_500rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_500rpart) <- c("5_fold_CV","DILI5")

dili5_1000rpart_mcc <- as.data.frame(rpart1000.5$pred)
dili5_1000rpart_mcc[,1] <- as.character(dili5_1000rpart_mcc[,1])
dili5_1000rpart_mcc[,2] <- as.character(dili5_1000rpart_mcc[,2])
dili5_1000rpart_mcc[dili5_1000rpart_mcc=="x1"] <- 1
dili5_1000rpart_mcc[dili5_1000rpart_mcc=="x0"] <- 0

dili5_1000rpart<-dili5_1000rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_1000rpart) <- c("5_fold_CV","DILI5")

dili6_100rpart_mcc <- as.data.frame(rpart100.6$pred)
dili6_100rpart_mcc[,1] <- as.character(dili6_100rpart_mcc[,1])
dili6_100rpart_mcc[,2] <- as.character(dili6_100rpart_mcc[,2])
dili6_100rpart_mcc[dili6_100rpart_mcc=="x1"] <- 1
dili6_100rpart_mcc[dili6_100rpart_mcc=="x0"] <- 0

dili6_100rpart<-dili6_100rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_100rpart) <- c("5_fold_CV","DILI6")

dili6_250rpart_mcc <- as.data.frame(rpart250.6$pred)
dili6_250rpart_mcc[,1] <- as.character(dili6_250rpart_mcc[,1])
dili6_250rpart_mcc[,2] <- as.character(dili6_250rpart_mcc[,2])
dili6_250rpart_mcc[dili6_250rpart_mcc=="x1"] <- 1
dili6_250rpart_mcc[dili6_250rpart_mcc=="x0"] <- 0

dili6_250rpart<-dili6_250rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_250rpart) <- c("5_fold_CV","DILI6")

dili6_500rpart_mcc <- as.data.frame(rpart500.6$pred)
dili6_500rpart_mcc[,1] <- as.character(dili6_500rpart_mcc[,1])
dili6_500rpart_mcc[,2] <- as.character(dili6_500rpart_mcc[,2])
dili6_500rpart_mcc[dili6_500rpart_mcc=="x1"] <- 1
dili6_500rpart_mcc[dili6_500rpart_mcc=="x0"] <- 0

dili6_500rpart<-dili6_500rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_500rpart) <- c("5_fold_CV","DILI6")

dili6_1000rpart_mcc <- as.data.frame(rpart1000.6$pred)
dili6_1000rpart_mcc[,1] <- as.character(dili6_1000rpart_mcc[,1])
dili6_1000rpart_mcc[,2] <- as.character(dili6_1000rpart_mcc[,2])
dili6_1000rpart_mcc[dili6_1000rpart_mcc=="x1"] <- 1
dili6_1000rpart_mcc[dili6_1000rpart_mcc=="x0"] <- 0

dili6_1000rpart<-dili6_1000rpart_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_1000rpart) <- c("5_fold_CV","DILI6")

dili1_100nb_mcc <- as.data.frame(NB100.1$pred)
dili1_100nb_mcc[,1] <- as.character(dili1_100nb_mcc[,1])
dili1_100nb_mcc[,2] <- as.character(dili1_100nb_mcc[,2])
dili1_100nb_mcc[dili1_100nb_mcc=="x1"] <- 1
dili1_100nb_mcc[dili1_100nb_mcc=="x0"] <- 0

dili1_100nb<-dili1_100nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_100nb) <- c("5_fold_CV","DILI1")

dili1_250nb_mcc <- as.data.frame(NB250.1$pred)
dili1_250nb_mcc[,1] <- as.character(dili1_250nb_mcc[,1])
dili1_250nb_mcc[,2] <- as.character(dili1_250nb_mcc[,2])
dili1_250nb_mcc[dili1_250nb_mcc=="x1"] <- 1
dili1_250nb_mcc[dili1_250nb_mcc=="x0"] <- 0

dili1_250nb<-dili1_250nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_250nb) <- c("5_fold_CV","DILI1")

dili1_500nb_mcc <- as.data.frame(NB500.1$pred)
dili1_500nb_mcc[,1] <- as.character(dili1_500nb_mcc[,1])
dili1_500nb_mcc[,2] <- as.character(dili1_500nb_mcc[,2])
dili1_500nb_mcc[dili1_500nb_mcc=="x1"] <- 1
dili1_500nb_mcc[dili1_500nb_mcc=="x0"] <- 0

dili1_500nb<-dili1_500nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_500nb) <- c("5_fold_CV","DILI1")

dili1_1000nb_mcc <- as.data.frame(NB1000.1$pred)
dili1_1000nb_mcc[,1] <- as.character(dili1_1000nb_mcc[,1])
dili1_1000nb_mcc[,2] <- as.character(dili1_1000nb_mcc[,2])
dili1_1000nb_mcc[dili1_1000nb_mcc=="x1"] <- 1
dili1_1000nb_mcc[dili1_1000nb_mcc=="x0"] <- 0

dili1_1000nb<-dili1_1000nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili1_1000nb) <- c("5_fold_CV","DILI1")

dili3_100nb_mcc <- as.data.frame(NB100.3$pred)
dili3_100nb_mcc[,1] <- as.character(dili3_100nb_mcc[,1])
dili3_100nb_mcc[,2] <- as.character(dili3_100nb_mcc[,2])
dili3_100nb_mcc[dili3_100nb_mcc=="x1"] <- 1
dili3_100nb_mcc[dili3_100nb_mcc=="x0"] <- 0

dili3_100nb<-dili3_100nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_100nb) <- c("5_fold_CV","DILI3")

dili3_250nb_mcc <- as.data.frame(NB250.3$pred)
dili3_250nb_mcc[,1] <- as.character(dili3_250nb_mcc[,1])
dili3_250nb_mcc[,2] <- as.character(dili3_250nb_mcc[,2])
dili3_250nb_mcc[dili3_250nb_mcc=="x1"] <- 1
dili3_250nb_mcc[dili3_250nb_mcc=="x0"] <- 0

dili3_250nb<-dili3_250nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_250nb) <- c("5_fold_CV","DILI3")

dili3_500nb_mcc <- as.data.frame(NB500.3$pred)
dili3_500nb_mcc[,1] <- as.character(dili3_500nb_mcc[,1])
dili3_500nb_mcc[,2] <- as.character(dili3_500nb_mcc[,2])
dili3_500nb_mcc[dili3_500nb_mcc=="x1"] <- 1
dili3_500nb_mcc[dili3_500nb_mcc=="x0"] <- 0

dili3_500nb<-dili3_500nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_500nb) <- c("5_fold_CV","DILI3")

dili3_1000nb_mcc <- as.data.frame(NB1000.3$pred)
dili3_1000nb_mcc[,1] <- as.character(dili3_1000nb_mcc[,1])
dili3_1000nb_mcc[,2] <- as.character(dili3_1000nb_mcc[,2])
dili3_1000nb_mcc[dili3_1000nb_mcc=="x1"] <- 1
dili3_1000nb_mcc[dili3_1000nb_mcc=="x0"] <- 0

dili3_1000nb<-dili3_1000nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili3_1000nb) <- c("5_fold_CV","DILI3")

dili5_100nb_mcc <- as.data.frame(NB100.5$pred)
dili5_100nb_mcc[,1] <- as.character(dili5_100nb_mcc[,1])
dili5_100nb_mcc[,2] <- as.character(dili5_100nb_mcc[,2])
dili5_100nb_mcc[dili5_100nb_mcc=="x1"] <- 1
dili5_100nb_mcc[dili5_100nb_mcc=="x0"] <- 0

dili5_100nb<-dili5_100nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_100nb) <- c("5_fold_CV","DILI5")

dili5_250nb_mcc <- as.data.frame(NB250.5$pred)
dili5_250nb_mcc[,1] <- as.character(dili5_250nb_mcc[,1])
dili5_250nb_mcc[,2] <- as.character(dili5_250nb_mcc[,2])
dili5_250nb_mcc[dili5_250nb_mcc=="x1"] <- 1
dili5_250nb_mcc[dili5_250nb_mcc=="x0"] <- 0

dili5_250nb<-dili5_250nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_250nb) <- c("5_fold_CV","DILI5")

dili5_500nb_mcc <- as.data.frame(NB500.5$pred)
dili5_500nb_mcc[,1] <- as.character(dili5_500nb_mcc[,1])
dili5_500nb_mcc[,2] <- as.character(dili5_500nb_mcc[,2])
dili5_500nb_mcc[dili5_500nb_mcc=="x1"] <- 1
dili5_500nb_mcc[dili5_500nb_mcc=="x0"] <- 0

dili5_500nb<-dili5_500nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_500nb) <- c("5_fold_CV","DILI5")

dili5_1000nb_mcc <- as.data.frame(NB1000.5$pred)
dili5_1000nb_mcc[,1] <- as.character(dili5_1000nb_mcc[,1])
dili5_1000nb_mcc[,2] <- as.character(dili5_1000nb_mcc[,2])
dili5_1000nb_mcc[dili5_1000nb_mcc=="x1"] <- 1
dili5_1000nb_mcc[dili5_1000nb_mcc=="x0"] <- 0

dili5_1000nb<-dili5_1000nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili5_1000nb) <- c("5_fold_CV","DILI5")

dili6_100nb_mcc <- as.data.frame(NB100.6$pred)
dili6_100nb_mcc[,1] <- as.character(dili6_100nb_mcc[,1])
dili6_100nb_mcc[,2] <- as.character(dili6_100nb_mcc[,2])
dili6_100nb_mcc[dili6_100nb_mcc=="x1"] <- 1
dili6_100nb_mcc[dili6_100nb_mcc=="x0"] <- 0

dili6_100nb<-dili6_100nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_100nb) <- c("5_fold_CV","DILI6")

dili6_250nb_mcc <- as.data.frame(NB250.6$pred)
dili6_250nb_mcc[,1] <- as.character(dili6_250nb_mcc[,1])
dili6_250nb_mcc[,2] <- as.character(dili6_250nb_mcc[,2])
dili6_250nb_mcc[dili6_250nb_mcc=="x1"] <- 1
dili6_250nb_mcc[dili6_250nb_mcc=="x0"] <- 0

dili6_250nb<-dili6_250nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_250nb) <- c("5_fold_CV","DILI6")

dili6_500nb_mcc <- as.data.frame(NB500.6$pred)
dili6_500nb_mcc[,1] <- as.character(dili6_500nb_mcc[,1])
dili6_500nb_mcc[,2] <- as.character(dili6_500nb_mcc[,2])
dili6_500nb_mcc[dili6_500nb_mcc=="x1"] <- 1
dili6_500nb_mcc[dili6_500nb_mcc=="x0"] <- 0

dili6_500nb<-dili6_500nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_500nb) <- c("5_fold_CV","DILI6")

dili6_1000nb_mcc <- as.data.frame(NB1000.6$pred)
dili6_1000nb_mcc[,1] <- as.character(dili6_1000nb_mcc[,1])
dili6_1000nb_mcc[,2] <- as.character(dili6_1000nb_mcc[,2])
dili6_1000nb_mcc[dili6_1000nb_mcc=="x1"] <- 1
dili6_1000nb_mcc[dili6_1000nb_mcc=="x0"] <- 0

dili6_1000nb<-dili6_1000nb_mcc%>%separate(Resample,c("cv","rep"),sep="\\.")%>%group_by(rep)%>%
    summarise(mcc=mccr(obs,pred))
colnames(dili6_1000nb) <- c("5_fold_CV","DILI6")

mccsvm100 <- merge(dili1_100svm, dili3_100svm, by="5_fold_CV")
mccsvm100 <- merge(mccsvm100, dili5_100svm, by="5_fold_CV")
mccsvm100 <- merge(mccsvm100, dili6_100svm, by="5_fold_CV")
colnames(mccsvm100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccsvm100, "../output_files/merge100svm_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

mccsvm250 <- merge(dili1_250svm, dili3_250svm, by="5_fold_CV")
mccsvm250 <- merge(mccsvm250, dili5_100svm, by="5_fold_CV")
mccsvm250 <- merge(mccsvm250, dili6_100svm, by="5_fold_CV")
colnames(mccsvm250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccsvm250, "../output_files/merge250svm_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccsvm500 <- merge(dili1_500svm, dili3_500svm, by="5_fold_CV")
mccsvm500 <- merge(mccsvm500, dili5_500svm, by="5_fold_CV")
mccsvm500 <- merge(mccsvm500, dili6_500svm, by="5_fold_CV")
colnames(mccsvm500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccsvm500, "../output_files/merge500svm_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccsvm1000 <- merge(dili1_1000svm, dili3_1000svm, by="5_fold_CV")
mccsvm1000 <- merge(mccsvm1000, dili5_1000svm, by="5_fold_CV")
mccsvm1000 <- merge(mccsvm1000, dili6_1000svm, by="5_fold_CV")
colnames(mccsvm1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccsvm1000, "../output_files/merge1000svm_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccrf100 <- merge(dili1_100rf, dili3_100rf, by="5_fold_CV")
mccrf100 <- merge(mccrf100, dili5_100rf, by="5_fold_CV")
mccrf100 <- merge(mccrf100, dili6_100rf, by="5_fold_CV")
colnames(mccrf100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccrf100, "../output_files/merge100rf_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

mccrf250 <- merge(dili1_250rf, dili3_250rf, by="5_fold_CV")
mccrf250 <- merge(mccrf250, dili5_100rf, by="5_fold_CV")
mccrf250 <- merge(mccrf250, dili6_100rf, by="5_fold_CV")
colnames(mccrf250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccrf250, "../output_files/merge250rf_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccrf500 <- merge(dili1_500rf, dili3_500rf, by="5_fold_CV")
mccrf500 <- merge(mccrf500, dili5_500rf, by="5_fold_CV")
mccrf500 <- merge(mccrf500, dili6_500rf, by="5_fold_CV")
colnames(mccrf500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccrf500, "../output_files/merge500rf_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccrf1000 <- merge(dili1_1000rf, dili3_1000rf, by="5_fold_CV")
mccrf1000 <- merge(mccrf1000, dili5_1000rf, by="5_fold_CV")
mccrf1000 <- merge(mccrf1000, dili6_1000rf, by="5_fold_CV")
colnames(mccrf1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccrf1000, "../output_files/merge1000rf_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccglm100 <- merge(dili1_100glm, dili3_100glm, by="5_fold_CV")
mccglm100 <- merge(mccglm100, dili5_100glm, by="5_fold_CV")
mccglm100 <- merge(mccglm100, dili6_100glm, by="5_fold_CV")
colnames(mccglm100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccglm100, "../output_files/merge100glm_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

mccglm250 <- merge(dili1_250glm, dili3_250glm, by="5_fold_CV")
mccglm250 <- merge(mccglm250, dili5_100glm, by="5_fold_CV")
mccglm250 <- merge(mccglm250, dili6_100glm, by="5_fold_CV")
colnames(mccglm250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccglm250, "../output_files/merge250glm_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccglm500 <- merge(dili1_500glm, dili3_500glm, by="5_fold_CV")
mccglm500 <- merge(mccglm500, dili5_500glm, by="5_fold_CV")
mccglm500 <- merge(mccglm500, dili6_500glm, by="5_fold_CV")
colnames(mccglm500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccglm500, "../output_files/merge500glm_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccglm1000 <- merge(dili1_1000glm, dili3_1000glm, by="5_fold_CV")
mccglm1000 <- merge(mccglm1000, dili5_1000glm, by="5_fold_CV")
mccglm1000 <- merge(mccglm1000, dili6_1000glm, by="5_fold_CV")
colnames(mccglm1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccglm1000, "../output_files/merge1000glm_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccrpart100 <- merge(dili1_100rpart, dili3_100rpart, by="5_fold_CV")
mccrpart100 <- merge(mccrpart100, dili5_100rpart, by="5_fold_CV")
mccrpart100 <- merge(mccrpart100, dili6_100rpart, by="5_fold_CV")
colnames(mccrpart100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccrpart100, "../output_files/merge100rpart_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

mccrpart250 <- merge(dili1_250rpart, dili3_250rpart, by="5_fold_CV")
mccrpart250 <- merge(mccrpart250, dili5_100rpart, by="5_fold_CV")
mccrpart250 <- merge(mccrpart250, dili6_100rpart, by="5_fold_CV")
colnames(mccrpart250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccrpart250, "../output_files/merge250rpart_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccrpart500 <- merge(dili1_500rpart, dili3_500rpart, by="5_fold_CV")
mccrpart500 <- merge(mccrpart500, dili5_500rpart, by="5_fold_CV")
mccrpart500 <- merge(mccrpart500, dili6_500rpart, by="5_fold_CV")
colnames(mccrpart500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccrpart500, "../output_files/merge500rpart_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccrpart1000 <- merge(dili1_1000rpart, dili3_1000rpart, by="5_fold_CV")
mccrpart1000 <- merge(mccrpart1000, dili5_1000rpart, by="5_fold_CV")
mccrpart1000 <- merge(mccrpart1000, dili6_1000rpart, by="5_fold_CV")
colnames(mccrpart1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccrpart1000, "../output_files/merge1000rpart_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccnb100 <- merge(dili1_100nb, dili3_100nb, by="5_fold_CV")
mccnb100 <- merge(mccnb100, dili5_100nb, by="5_fold_CV")
mccnb100 <- merge(mccnb100, dili6_100nb, by="5_fold_CV")
colnames(mccnb100) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccnb100, "../output_files/merge100nb_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv", quote=FALSE, row.names=FALSE)

mccnb250 <- merge(dili1_250nb, dili3_250nb, by="5_fold_CV")
mccnb250 <- merge(mccnb250, dili5_100nb, by="5_fold_CV")
mccnb250 <- merge(mccnb250, dili6_100nb, by="5_fold_CV")
colnames(mccnb250) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccnb250, "../output_files/merge250nb_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccnb500 <- merge(dili1_500nb, dili3_500nb, by="5_fold_CV")
mccnb500 <- merge(mccnb500, dili5_500nb, by="5_fold_CV")
mccnb500 <- merge(mccnb500, dili6_500nb, by="5_fold_CV")
colnames(mccnb500) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccnb500, "../output_files/merge500nb_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)

mccnb1000 <- merge(dili1_1000nb, dili3_1000nb, by="5_fold_CV")
mccnb1000 <- merge(mccnb1000, dili5_1000nb, by="5_fold_CV")
mccnb1000 <- merge(mccnb1000, dili6_1000nb, by="5_fold_CV")
colnames(mccnb1000) <- c("CAM_ID","DILI1","DILI3","DILI5","DILI6")
write.csv(mccnb1000, "../output_files/merge1000nb_p3p4p5p6p7p8-crossvalidation-camda2020-UND_MAR2021.csv.csv", quote=FALSE, row.names=FALSE)


