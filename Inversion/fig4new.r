#!/usr/bin/R

library(ggpubr)
library(gtools)
library(dplyr)
library(reshape2)
library(SeqVarTools)

# file: "lg18.ld","lg18_inv.vcf","lg21.ld","lg21_inv.vcf","lg26.ld","lg26_inv.vcf"
hplots <- list()
splots <- list()
zplots <- list()
lplots <- list()
######################################################################
#palette =c("#ff0000","#ff9900","#00ff00","#0000ff","#00ffff","#9900ff"),
#ch,hl,hz,lk,th,yl
######################################################################
#getleg
abbr <- c("CH", "HL", "HZ", "LK", "TH", "YL")
names(abbr) <- c("CH", "HL", "HZ", "LK", "TH", "YL")
vcf <- "26.temp"
cat("  Processing vcf", vcf, "...")
seqVCF2GDS(vcf, ".tmp.gds", parallel=8,verbose=FALSE)
gds <- seqOpen(".tmp.gds")
pc  <- pca(gds, eigen.cnt=length(seqGetData(gds, "sample.id")))
pc.va <- pc$eigenval
set.seed(222) # for reproducible
# partinate to 3 clusters by kmeans.
clst <- kmeans(pc$eigenvect[,1], 3, iter.max=100, nstart=100)
clst.o  <- order(clst$centers)
hap1 <- names(clst$cluster[which(clst$cluster==clst.o[1])])
hap2 <- names(clst$cluster[which(clst$cluster==clst.o[2])])
hap3 <- names(clst$cluster[which(clst$cluster==clst.o[3])])
write.table(hap1, file=paste0(vcf,".hap1"), quote=F, row.names=F, col.names=F)
write.table(hap2, file=paste0(vcf,".hap2"), quote=F, row.names=F, col.names=F)
write.table(hap3, file=paste0(vcf,".hap3"), quote=F, row.names=F, col.names=F)
# calc het.
het <- list(
"AA"=heterozygosity(gds, margin="by.sample", use.names=T)[hap1], 
"AB"=heterozygosity(gds, margin="by.sample", use.names=T)[hap2], 
"BB"=heterozygosity(gds, margin="by.sample", use.names=T)[hap3])
md <- melt(het, varnames=c("het", "hap"))
names(md) <- c("het","hap")

# boxplot.
boxp <- ggboxplot(md, x="hap", y="het", 
color = "hap",add = "jitter", shape="hap", 
xlab="", ylab="Heterozygosity",title="",palette = "lancet") +
theme(legend.position="none", plot.title=element_text(hjust=0.5))

# pca dataframe.
v <- rep(c("AA","AB","BB"), c(length(hap1),length(hap2),length(hap3)))
names(v) <- c(hap1,hap2,hap3)
dat <- data.frame(
	"ind"=row.names(pc$eigenvect), 
	"clusters"=v[row.names(pc$eigenvect)],
	"Populations"=abbr[substr(row.names(pc$eigenvect),1,2)],
	"PC1"=pc$eigenvect[,1], "PC2"=pc$eigenvect[,2],
	stringsAsFactors = FALSE)
dat$Populations <- factor(dat$Populations, 
	levels = c("CH", "HL", "HZ", "LK", "TH", "YL"))

# pca scatter plot.
lab1<-paste("PC1"," (",sprintf("%.2f%%", 100*pc.va[1]/sum(pc.va)),")", sep="")
lab2<-paste("PC2"," (",sprintf("%.2f%%", 100*pc.va[2]/sum(pc.va)),")", sep="")
scap<-ggscatter(dat,x="PC1",y="PC2",size=2, font.label = c(0.1, "plain"),
    color="Populations",shape="clusters",
    title="",xlab=lab1,ylab=lab2,palette =c("#ff0000","#ff9900","#00ff00","#0000ff","#00ffff","#9900ff")) +
    theme(
    legend.title=element_blank(),
	legend.position="left",
    plot.title=element_text(hjust=0.5))
	
seqClose(gds)
cat("done.\n")
leg <- get_legend(scap)
dev.off()











######################################################################
# lg18
## ld
ld <- read.table("lg18",header=T)
ld$POS1 <- ld$POS1*50000/1000000
ld$POS2 <- ld$POS2*50000/1000000
ldp1 <- ggplot(ld,aes(x=POS1,y=POS2)) +
theme_bw() +
geom_tile(aes(fill=LD)) + scale_fill_continuous(type = "viridis" )+
labs(title="LG18",x="Chromosome 18 (Mbp)",y="Chromosome 18 (Mbp)") +
theme(plot.title = element_text(size=30,hjust = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


abbr <- c("CH", "HZ", "LK", "TH", "YL")
names(abbr) <- c("CH", "HZ", "LK", "TH", "YL")
vcf <- "lg18_inv.vcf"
cat("  Processing vcf", vcf, "...")
seqVCF2GDS(vcf, ".tmp.gds", parallel=8,verbose=FALSE)
gds <- seqOpen(".tmp.gds")
pc  <- pca(gds, eigen.cnt=length(seqGetData(gds, "sample.id")))
pc.va <- pc$eigenval
set.seed(222) # for reproducible
# partinate to 3 clusters by kmeans.
clst <- kmeans(pc$eigenvect[,1], 3, iter.max=100, nstart=100)
clst.o  <- order(clst$centers)
hap1 <- names(clst$cluster[which(clst$cluster==clst.o[1])])
hap2 <- names(clst$cluster[which(clst$cluster==clst.o[2])])
hap3 <- names(clst$cluster[which(clst$cluster==clst.o[3])])
write.table(hap1, file=paste0(vcf,".hap1"), quote=F, row.names=F, col.names=F)
write.table(hap2, file=paste0(vcf,".hap2"), quote=F, row.names=F, col.names=F)
write.table(hap3, file=paste0(vcf,".hap3"), quote=F, row.names=F, col.names=F)
# calc het.
het <- list(
"AA"=heterozygosity(gds, margin="by.sample", use.names=T)[hap1], 
"AB"=heterozygosity(gds, margin="by.sample", use.names=T)[hap2], 
"BB"=heterozygosity(gds, margin="by.sample", use.names=T)[hap3])
md <- melt(het, varnames=c("het", "hap"))
names(md) <- c("het","hap")

# boxplot.
boxp1 <- ggboxplot(md, x="hap", y="het", 
color = "hap",add = "jitter", shape="hap", 
xlab="", ylab="Heterozygosity",title="",palette = "lancet") +
theme(legend.position="none", plot.title=element_text(hjust=0.5))

# pca dataframe.
v <- rep(c("AA","AB","BB"), c(length(hap1),length(hap2),length(hap3)))
names(v) <- c(hap1,hap2,hap3)
dat <- data.frame(
	"ind"=row.names(pc$eigenvect), 
	"clusters"=v[row.names(pc$eigenvect)],
	"Populations"=abbr[substr(row.names(pc$eigenvect),1,2)],
	"PC1"=pc$eigenvect[,1], "PC2"=pc$eigenvect[,2],
	stringsAsFactors = FALSE)
dat$Populations <- factor(dat$Populations, 
	levels = c("CH", "HZ", "LK", "TH", "YL"))

# pca scatter plot.
lab1<-paste("PC1"," (",sprintf("%.2f%%", 100*pc.va[1]/sum(pc.va)),")", sep="")
lab2<-paste("PC2"," (",sprintf("%.2f%%", 100*pc.va[2]/sum(pc.va)),")", sep="")
scap1<-ggscatter(dat,x="PC1",y="PC2",size=2, font.label = c(0.1, "plain"),
    color="Populations",shape="clusters",
    title="",xlab=lab1,ylab=lab2,palette =c("#ff0000","#00ff00","#0000ff","#00ffff","#9900ff")) +
    theme(
    legend.title=element_blank(),
	legend.position="left",
    plot.title=element_text(hjust=0.5))
	
seqClose(gds)
cat("done.\n")

#barplot
aa <- data.frame(ind=hap1,clu="AA")
ab <- data.frame(ind=hap2,clu="AB")
bb <- data.frame(ind=hap3,clu="BB")
h_d <- h_d <- rbind(aa,ab,bb)
h_d <- mutate(h_d,pop=substring(ind,1,2))
hdcount <- count(h_d,pop,clu)
hdallct <- hdcount %>%group_by(pop)%>% mutate(fre= n/sum(n))
barp1 <- ggbarplot(hdallct,x="pop",y="fre",xlab=FALSE,ylab="Frequency",legend.title="",fill="clu",color="clu",palette = "lancet",order=c("LK","YL","CH","HZ","TH"))
######################################################################
# lg21
## ld
ld <- read.table("lg21",header=T)
ld$POS1 <- ld$POS1*50000/1000000
ld$POS2 <- ld$POS2*50000/1000000
ldp2 <- ggplot(ld,aes(x=POS1,y=POS2)) +
theme_bw() +
geom_tile(aes(fill=LD)) + scale_fill_continuous(type = "viridis" )+
labs(title="LG21",x="Chromosome 21 (Mbp)",y="Chromosome 21 (Mbp)") +
theme(plot.title = element_text(size=30,hjust = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


abbr <- c("CH", "HL", "HZ", "TH")
names(abbr) <- c("CH", "HL", "HZ", "TH")
vcf <- "lg21_inv.vcf"
cat("  Processing vcf", vcf, "...")
seqVCF2GDS(vcf, ".tmp.gds", parallel=8,verbose=FALSE)
gds <- seqOpen(".tmp.gds")
pc  <- pca(gds, eigen.cnt=length(seqGetData(gds, "sample.id")))
pc.va <- pc$eigenval
set.seed(222) # for reproducible
# partinate to 3 clusters by kmeans.
clst <- kmeans(pc$eigenvect[,1], 3, iter.max=100, nstart=100)
clst.o  <- order(clst$centers)
hap1 <- names(clst$cluster[which(clst$cluster==clst.o[1])])
hap2 <- names(clst$cluster[which(clst$cluster==clst.o[2])])
hap3 <- names(clst$cluster[which(clst$cluster==clst.o[3])])
write.table(hap1, file=paste0(vcf,".hap1"), quote=F, row.names=F, col.names=F)
write.table(hap2, file=paste0(vcf,".hap2"), quote=F, row.names=F, col.names=F)
write.table(hap3, file=paste0(vcf,".hap3"), quote=F, row.names=F, col.names=F)
# calc het.
het <- list(
"AA"=heterozygosity(gds, margin="by.sample", use.names=T)[hap1], 
"AB"=heterozygosity(gds, margin="by.sample", use.names=T)[hap2], 
"BB"=heterozygosity(gds, margin="by.sample", use.names=T)[hap3])
md <- melt(het, varnames=c("het", "hap"))
names(md) <- c("het","hap")

# boxplot.
boxp2 <- ggboxplot(md, x="hap", y="het", 
color = "hap",add = "jitter", shape="hap", 
xlab="",ylab="", title="",palette = "lancet") +
theme(legend.position="none", plot.title=element_text(hjust=0.5))

# pca dataframe.
v <- rep(c("AA","AB","BB"), c(length(hap1),length(hap2),length(hap3)))
names(v) <- c(hap1,hap2,hap3)
dat <- data.frame(
	"ind"=row.names(pc$eigenvect), 
	"clusters"=v[row.names(pc$eigenvect)],
	"Populations"=abbr[substr(row.names(pc$eigenvect),1,2)],
	"PC1"=pc$eigenvect[,1], "PC2"=pc$eigenvect[,2],
	stringsAsFactors = FALSE)
dat$Populations <- factor(dat$Populations, 
	levels = c("CH", "HL", "HZ", "TH"))

# pca scatter plot.
lab1<-paste("PC1"," (",sprintf("%.2f%%", 100*pc.va[1]/sum(pc.va)),")", sep="")
lab2<-paste("PC2"," (",sprintf("%.2f%%", 100*pc.va[2]/sum(pc.va)),")", sep="")
scap2<-ggscatter(dat,x="PC1",y="PC2",size=2, font.label = c(0.1, "plain"),
    color="Populations",shape="clusters",
    title="",xlab=lab1,ylab=lab2,palette =c("#ff0000","#ff9900","#00ff00","#00ffff")) +
    theme(
    legend.title=element_blank(),
	legend.position="left",
    plot.title=element_text(hjust=0.5))
	
seqClose(gds)
cat("done.\n")
#barplot
aa <- data.frame(ind=hap1,clu="AA")
ab <- data.frame(ind=hap2,clu="AB")
bb <- data.frame(ind=hap3,clu="BB")
h_d <- h_d <- rbind(aa,ab,bb)
h_d <- mutate(h_d,pop=substring(ind,1,2))
hdcount <- count(h_d,pop,clu)
hdallct <- hdcount %>%group_by(pop)%>% mutate(fre= n/sum(n))
barp2 <- ggbarplot(hdallct,x="pop",y="fre",xlab=FALSE,ylab=FALSE,legend.title="",fill="clu",color="clu",palette = "lancet",order=c("HL","CH","HZ","TH"))
######################################################################
# lg26
## ld
ld <- read.table("lg26",header=T)
ld$POS1 <- ld$POS1*50000/1000000
ld$POS2 <- ld$POS2*50000/1000000
ldp3 <- ggplot(ld,aes(x=POS1,y=POS2)) +
theme_bw() +
geom_tile(aes(fill=LD)) + scale_fill_continuous(type = "viridis" )+
labs(title="LG26",x="Chromosome 26 (Mbp)",y="Chromosome 26 (Mbp)") +
theme(plot.title = element_text(size=30,hjust = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


abbr <- c("CH", "HL", "HZ", "TH")
names(abbr) <- c("CH", "HL", "HZ", "TH")
vcf <- "lg26_inv.vcf"
cat("  Processing vcf", vcf, "...")
seqVCF2GDS(vcf, ".tmp.gds", parallel=8,verbose=FALSE)
gds <- seqOpen(".tmp.gds")
pc  <- pca(gds, eigen.cnt=length(seqGetData(gds, "sample.id")))
pc.va <- pc$eigenval
set.seed(222) # for reproducible
# partinate to 3 clusters by kmeans.
clst <- kmeans(pc$eigenvect[,1], 3, iter.max=100, nstart=100)
clst.o  <- order(clst$centers)
hap1 <- names(clst$cluster[which(clst$cluster==clst.o[1])])
hap2 <- names(clst$cluster[which(clst$cluster==clst.o[2])])
hap3 <- names(clst$cluster[which(clst$cluster==clst.o[3])])
write.table(hap1, file=paste0(vcf,".hap1"), quote=F, row.names=F, col.names=F)
write.table(hap2, file=paste0(vcf,".hap2"), quote=F, row.names=F, col.names=F)
write.table(hap3, file=paste0(vcf,".hap3"), quote=F, row.names=F, col.names=F)
# calc het.
het <- list(
"AA"=heterozygosity(gds, margin="by.sample", use.names=T)[hap1], 
"AB"=heterozygosity(gds, margin="by.sample", use.names=T)[hap2], 
"BB"=heterozygosity(gds, margin="by.sample", use.names=T)[hap3])
md <- melt(het, varnames=c("het", "hap"))
names(md) <- c("het","hap")

# boxplot.
boxp3 <- ggboxplot(md, x="hap", y="het", 
color = "hap",add = "jitter", shape="hap", 
xlab="", ylab="",title="",palette = "lancet") +
theme(legend.position="none", plot.title=element_text(hjust=0.5))

# pca dataframe.
v <- rep(c("AA","AB","BB"), c(length(hap1),length(hap2),length(hap3)))
names(v) <- c(hap1,hap2,hap3)
dat <- data.frame(
	"ind"=row.names(pc$eigenvect), 
	"clusters"=v[row.names(pc$eigenvect)],
	"Populations"=abbr[substr(row.names(pc$eigenvect),1,2)],
	"PC1"=pc$eigenvect[,1], "PC2"=pc$eigenvect[,2],
	stringsAsFactors = FALSE)
dat$Populations <- factor(dat$Populations, 
	levels = c("CH", "HL", "HZ", "TH"))

# pca scatter plot.
lab1<-paste("PC1"," (",sprintf("%.2f%%", 100*pc.va[1]/sum(pc.va)),")", sep="")
lab2<-paste("PC2"," (",sprintf("%.2f%%", 100*pc.va[2]/sum(pc.va)),")", sep="")
scap3<-ggscatter(dat,x="PC1",y="PC2",size=2, font.label = c(0.1, "plain"),
    color="Populations",shape="clusters",
    title="",xlab=lab1,ylab=lab2,palette =c("#ff0000","#ff9900","#00ff00","#00ffff")) +
    theme(
    legend.title=element_blank(),
	legend.position="left",
    plot.title=element_text(hjust=0.5))
	
seqClose(gds)
cat("done.\n")
#barplot
aa <- data.frame(ind=hap1,clu="AA")
ab <- data.frame(ind=hap2,clu="AB")
bb <- data.frame(ind=hap3,clu="BB")
h_d <- h_d <- rbind(aa,ab,bb)
h_d <- mutate(h_d,pop=substring(ind,1,2))
hdcount <- count(h_d,pop,clu)
hdallct <- hdcount %>%group_by(pop)%>% mutate(fre= n/sum(n))
barp3 <- ggbarplot(hdallct,x="pop",y="fre",xlab=FALSE,ylab=FALSE,legend.title="",fill="clu",color="clu",palette = "lancet",order=c("HL","CH","HZ","TH"))
######################################################################

pdf("fig4new.pdf", width=14, height=17,onefile=F)
ldp <- ggarrange(ldp1,ldp2,ldp3, ncol = 3, nrow=1, common.legend = T, legend="right") 
lds <- ggarrange(scap1,scap2,scap3, ncol = 3, nrow=1, common.legend = T,legend.grob = leg, legend="right")

ldb <- ggarrange(boxp1,boxp2,boxp3, ncol = 3, nrow=1)
lda <- ggarrange(barp1,barp2,barp3,ncol = 3, nrow=1, common.legend = T, legend="right")
#ggarrange(ldp,lds,ldb,lda, labels =c("a","b","c","d"),hjust = -1,font.label = list(size =32, face = "bold", color ="black"),ncol = 1, nrow=4)
 ggarrange(ldp,lds,ldb,lda, labels ="AUTO",ncol = 1, nrow=4)
dev.off()

