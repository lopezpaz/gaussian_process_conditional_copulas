library(gplots)
library(tikzDevice)

tikz("bars.tex", standAlone=TRUE,height=1.5, width=8)

par(mfrow=c(1,7))
par(mar=c(2.5,2,1,0.5))
i <- 1
for(d in c("cloud", "glass", "housing", "jura", "shuttle", "weather", "stocks")) {
x <- read.table("results")
y <- subset(x, V1==d)

if(i == 1) {
  plotCI(x=(1:nrow(y))+1,y=y[,2],uiw=y[,3], 
         ylim=c(min(y[,2],y[,4])-max(y[,3],y[,5]), 
                max(y[,2],y[,4])+max(y[,3],y[,5])),
         ylab="Log-Likelihood", xlab="", main=d, lwd=2)
} else {
  plotCI(x=(1:nrow(y))+1,y=y[,2],uiw=y[,3], 
         ylim=c(min(y[,2],y[,4])-max(y[,3],y[,5]), 
                max(y[,2],y[,4])+max(y[,3],y[,5])),
         ylab="", xlab="", main=d, lwd=2)
}

plotCI(x=(1:nrow(y))+1,y=y[,4],uiw=y[,5],add=TRUE,col="red", lwd=2,pch=2)
i <- i + 1
}
  dev.off(); system("pdflatex bars.tex")
