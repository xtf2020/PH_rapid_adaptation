library(pcadapt)
library(qvalue)
filename <- read.pcadapt('input.bed', type = "bed")
res <- pcadapt(filename, K = 20, LD.clumping = NULL)
pdf("k20.pdf")
plot(res, option = "screeplot")
dev.off()
k<-2
res <- pcadapt(filename, K = as.numeric(k), LD.clumping = NULL)
qval <- qvalue(res$pvalues)$qvalues
outliers <- which(qval <= 0.1)

write.table(outliers, "pcadapt.outliers.id", col.names=F, row.names=F, quote=F)
