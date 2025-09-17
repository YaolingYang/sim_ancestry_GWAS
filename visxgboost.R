xgbpredres=read.csv("xgbpred_testdata.csv",row.names=1)

birthplaces=read.table("birthplace_UKB35k.txt",header=T)
rownames(birthplaces)=birthplaces[,1]
birthplaces=birthplaces[rownames(xgbpredres),]

cols=as.numeric(as.factor(birthplaces$birthcat))
xlim=range(xgbpredres[,6])
ylim=range(xgbpredres[,5])
ylim[1]=ylim[1]+30000
ylim[2]=ylim[2]-150000

## Plot for Figure 1
cex=1
cex.text=6
#pdf("xgbvis_k.pdf",height=12,width=16)
png("xgbvis_k.png",height=1200,width=2000)
par(mfrow=c(1,3),mar=c(1,4,1,4))
par(xpd=NA)
plot(xgbpredres[,2],xgbpredres[,1],pch=19,cex=cex,xlim=xlim,
     ylim=ylim,col=cols,main="",axes=F,xlab="",ylab="")
mtext("d",adj=0,cex=cex.text,line=-5)
text(600000,600000,"prediction\nusing\nHCs",cex=cex.text)
plot(xgbpredres[,4],xgbpredres[,3],pch=19,cex=cex,xlim=xlim,
     ylim=ylim,col=cols,main="",axes=F,xlab="",ylab="")
text(600000,600000,"prediction\nusing\nPCs",cex=cex.text)
plot(xgbpredres[,6],xgbpredres[,5],pch=19,cex=cex,
     col=cols,xlim=xlim,ylim=ylim,main="",axes=F,xlab="",ylab="")
text(600000,600000,"observed\nUK\nbirthplace",cex=cex.text)
dev.off()

## Completing the plot for figure 1

library(magick)
library(grid)
library(gridExtra)
## "Figure1.pdf" is the top part of the figure
g1 <- rasterGrob(as.raster(image_read_pdf("Figure1.pdf")))
g2 <- rasterGrob(as.raster(image_read("xgbvis_k.png")))
pdf("Figure1full.pdf", width = 15, height = 15)
grid.arrange(g1, g2, ncol = 1, heights = c(10,5))
dev.off()

## Just the HCs
cols=as.numeric(as.factor(birthplaces$birthcat))
tcols=col2rgb(cols,alpha = TRUE)
tcols["alpha",]=50
colsalpha=rgb(tcols[1,],tcols[2,],tcols[3,],tcols[4,],maxColorValue=255)

png("xgbvis_striking.png",height=1200*2,width=800*2)
plot(xgbpredres[,2],xgbpredres[,1],pch=19,cex=1.5,xlim=xlim,
     ylim=ylim,col=colsalpha,main="",axes=F,xlab="",ylab="")
dev.off()
