#!/usr/bin/R

library(ggpubr)
library(gtools)
library(dplyr)
library(reshape2)
library(SeqVarTools)

# pop abbr
abbr <- c("CH", "HL", "HZ", "LK", "TH", "YL")
names(abbr) <- c("CH", "HL", "HZ", "LK", "TH", "YL")

argv    <- commandArgs(TRUE)
vcf_dir <- argv[1]
file_names <- list.files(path=vcf_dir, pattern='.vcf$')
file_names <- mixedsort(file_names)
plots <- list()
plots1 <- list()

a=ceiling(length(file_names)/3)
b=ceiling(length(file_names)/1)

for (i in 1:length(file_names)) {
    vcf <- file_names[i]
	cat("  Processing vcf", vcf, "...")
    LG  <- gsub(".recode.vcf","",vcf)
	name<- strsplit(LG, "_")
    seqVCF2GDS(paste(vcf_dir,vcf,sep="/"), ".tmp.gds", parallel=8,verbose=FALSE)
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
	hapbox <- ggboxplot(md, x="hap", y="het", 
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
    p<-ggscatter(dat,x="PC1",y="PC2",size=2, font.label = c(0.1, "plain"),
    color="Populations",shape="clusters",
    title="",xlab=lab1,ylab=lab2,palette =c("#ff0000","#ff9900","#00ff00","#0000ff","#00ffff","#9900ff")) +
    theme(
    legend.title=element_blank(),
	legend.position="left",
    plot.title=element_text(hjust=0.5))
    #barplot
	aa <- data.frame(ind=hap1,clu="AA")
ab <- data.frame(ind=hap2,clu="AB")
bb <- data.frame(ind=hap3,clu="BB")
h_d <- h_d <- rbind(aa,ab,bb)
h_d <- mutate(h_d,pop=substring(ind,1,2))
hdcount <- count(h_d,pop,clu)
hdallct <- hdcount %>%group_by(pop)%>% mutate(fre= n/sum(n))
barp <- ggbarplot(hdallct,x="pop",y="fre",xlab=FALSE,ylab="Frequency",legend.title="",fill="clu",color="clu",palette = "lancet",order=c("HL","CH","HZ","TH","LK","YL"))
	
    if (i%%3 != 1) {p <- p + theme(legend.position="none")}
    plots1[[i]] <- p
	cmb <- ggarrange(p+labs(title="")+ theme(legend.position="none"), 
		hapbox+labs(title="")+ theme(legend.position="none"), 
		barp+labs(title="")+ theme(legend.position="none"),ncol = 3)
	cmb <- annotate_figure(cmb, 
		top = text_grob(paste0(name[[1]][1], 
		": Median.LD = ", name[[1]][4], 
		" No. Loci = ", name[[1]][5])))
	plots[[i]] <- cmb
	
	
    # clean files.
	seqClose(gds)
	cat("done.\n")
}

leg <- get_legend(plots1[[1]])
dev.off()

pdf("inv_all_pca_box_bar.pdf", width=14, height=b*3)
ggarrange(plotlist=plots, ncol = 1, nrow = b,
	labels="AUTO",
	common.legend = T, legend.grob = leg, legend="left") 
dev.off()
