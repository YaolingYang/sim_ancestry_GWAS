## This file performs xgboost prediction from either HCs or PCs as documented in
## [Yang, Y. & Lawson, D.J. From individuals to ancestries: towards attributing trait variation to haplotypes.
## PLoS Genet (2025) to appear. Available from medRxiv (2025). doi:10.1101/2024.03.13.24304206](https://www.medrxiv.org/content/10.1101/2025.03.13.25323895v1)

## Input files: "UKB_id_QC.txt": single column file consisting of the UK Biobank individual identifiers
## ethnic.txt: mapping ids from the above to self-reported ethnicity
## data_birthplaces.txt: UKB data extract in the format (with header) f.eid   f.54.0.0        f.129.0.0       f.130.0.0       f.20115.0.0     f.21000.0.0     f.22006.0.0
##      (where f.eid can be matched with UKB_id_QC.txt)

## "flashpca_200.pcs.gz": File as output by flashpca, of format (with header) FID     IID     PC1     PC2     PC3     PC4     PC5 ...
## "150HCs_UKB407k_log.csv.gz": Output of pbwt -PaintSparse, read as a matrix and then logged, of format (without header) IID HC1 HC2 ... 

##begin analysing UKB PCs
a=read.table("UKB_id_QC.txt")
b=read.table("ethnic.txt")
c=match(a[,1],b[,1])
matchdata=b[c,]
matchdata[which(is.na(matchdata[,2])), 2] <- "Other/unknown" #grey
matchdata[matchdata[,2] == -1, 2] <- "Other/unknown"
matchdata[matchdata[,2] == -3, 2] <- "Other/unknown" 
matchdata[matchdata[,2] == 1, 2] <- "Any other white background" #red
matchdata[matchdata[,2] == 1001, 2] <- "British" #red
matchdata[matchdata[,2] == 1002, 2] <- "Irish" #red
matchdata[matchdata[,2] == 1003, 2] <- "Any other white background" 
matchdata[matchdata[,2] == 2, 2] <- "Any other mixed background" #grey
matchdata[matchdata[,2] == 2001, 2] <- "White and Black Caribbean" #yellow
matchdata[matchdata[,2] == 2002, 2] <- "White and Black African" #yellow
matchdata[matchdata[,2] == 2003, 2] <- "White and Asian" #green
matchdata[matchdata[,2] == 2004, 2] <- "Any other mixed background"
matchdata[matchdata[,2] == 3, 2] <- "Any other Asian background" #blue
matchdata[matchdata[,2] == 3001, 2] <- "Indian" #blue
matchdata[matchdata[,2] == 3002, 2] <- "Pakistani" #blue
matchdata[matchdata[,2] == 3003, 2] <- "Bangladeshi" #blue
matchdata[matchdata[,2] == 3004, 2] <- "Any other Asian background"
matchdata[matchdata[,2] == 4, 2] <- "Any other Black background" #orange
matchdata[matchdata[,2] == 4001, 2] <- "Caribbean" #orange
matchdata[matchdata[,2] == 4002, 2] <- "African" #orange
matchdata[matchdata[,2] == 4003, 2] <- "Any other Black background"
matchdata[matchdata[,2] == 5, 2] <- "Chinese" #purple
matchdata[matchdata[,2] == 6, 2] <- "Other/unknown"

## READ THE PCS

GC=read.table("flashpca_200.pcs.gz",header = T)
c2=match(matchdata[,1],GC[,1])
GC=GC[c2,3:152]

## READ THE HCs
HC=read.table("150HCs_UKB407k_log.csv.gz",sep=",",header=F)
HC=HC[,-1]

HC$ethnicity <- matchdata[,2]
library(ggplot2)
ethnicity_colors <- c(
  "Other/unknown" = "#808080", 
  "Any other white background" = "#FF0000", 
  "British" = "#B22222", 
  "Irish" = "#8B0000", 
  "Any other mixed background" = "#808080", 
  "White and Black Caribbean" = "#FFFF00", 
  "White and Black African" = "#FFD700", 
  "White and Asian" = "orange",  
  "Any other Asian background" = "#0000FF", 
  "Indian" = "#1E90FF", 
  "Pakistani" = "#4169E1", 
  "Bangladeshi" = "#00008B", 
  "Any other Black background" = "#90EE90", 
  "Caribbean" = "#32CD32", 
  "African" = "#006400", 
  "Chinese" = "#800080"  
)
HC$ethnicity <- factor(HC$ethnicity, levels = names(ethnicity_colors))

## extract the individuals born in the uk
## So that we can use xgboost to fit birthplaces (east and north coordianates) against 150HCs/GCs
birth=read.table("data_birthplaces.txt",header=T)[,c(1,3,4,6)]
birth=birth[birth[,4]%in% c(1,1001,1002,1003),]
birth=birth[-which(is.na(birth[,1])),]
birth=birth[-which(is.na(birth[,2])),]
birth=birth[-which(birth[,2]==-1),]
c3=match(matchdata[,1],birth[,1])
uk=which(!is.na(c3))
#351213 individuals in the UK
birth=birth[c3[uk],-1]
HC_uk=HC[uk,1:150]
GC_uk=GC[uk,1:150]

### Perform xgboost
library(xgboost)
set.seed(1)
train_sample=sample(1:nrow(birth),ceiling(nrow(birth)*0.8))
data=cbind(birth[,1],HC_uk)
colnames(data)[1]="east"

xgb <- xgboost(data = as.matrix(data[train_sample, -1]),
               label = data$east[train_sample], nrounds = 45)
xgb_pred_east_HC <- predict(xgb, as.matrix(data[-train_sample, -1]))

data=cbind(birth[,2],HC_uk)
colnames(data)[1]="north"
xgb <- xgboost(data = as.matrix(data[train_sample, -1]),
               label = data$north[train_sample], nrounds = 45)
xgb_pred_north_HC <- predict(xgb, as.matrix(data[-train_sample, -1]))

data=cbind(birth[,1],GC_uk)
colnames(data)[1]="east"
xgb <- xgboost(data = as.matrix(data[train_sample, -1]),
               label = data$east[train_sample], nrounds = 45)
xgb_pred_east_GC <- predict(xgb, as.matrix(data[-train_sample, -1]))

data=cbind(birth[,2],GC_uk)
colnames(data)[1]="north"
xgb <- xgboost(data = as.matrix(data[train_sample, -1]),
               label = data$north[train_sample], nrounds = 45)
xgb_pred_north_GC <- predict(xgb, as.matrix(data[-train_sample, -1]))

## Tidy up and save
xgbpred=matrix(0,nrow=dim(data)[1] - length(train_sample),ncol=7)
colnames(xgbpred)=c("east_HC","north_HC","east_GC","north_GC","east","north","idx")
rownames(xgbpred)=a[uk,,drop=FALSE][-train_sample,1]
xgbpred[,1] = xgb_pred_east_HC
xgbpred[,2] = xgb_pred_north_HC
xgbpred[,3] = xgb_pred_east_GC
xgbpred[,4] = xgb_pred_north_GC
xgbpred[,5] = birth[-train_sample,1]
xgbpred[,6] = birth[-train_sample,2]
xgbpred[,7] = (1:dim(data)[1])[-train_sample]
write.csv(xgbpred,file="xgbpred_testdata.csv")
