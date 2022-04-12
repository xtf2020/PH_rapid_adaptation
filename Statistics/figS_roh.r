library(ggpubr)
library(dplyr)


d <- read.table("roh.indiv",header=T)
d <- mutate(d,p=substr(IID,1,2))
com <- list (c("HL","CH"),c("HL","HZ"),c("HL","TH"))
d$KB = d$KB/1000
n <- ggboxplot(d,x="p",y="NSEG",color="p",palette="npg",order=c("HL","CH","HZ","TH","LK","YL"),xlab="",ylab="number of ROHs",legend="none")+\
stat_compare_means(comparisons=com,method="wilcox.test")
l <- ggboxplot(d,x="p",y="KB",color="p",palette="npg",order=c("HL","CH","HZ","TH","LK","YL"),xlab="",ylab="length of ROHs (MB)",legend="none")+\
stat_compare_means(comparisons=com,method="wilcox.test")

jpeg("figS_roh.jpg",width=600,height=1200)
ggarrange(n,l,ncol=1)
dev.off()
