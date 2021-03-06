library(tidyr)
library(dplyr)
library(ggplot2)

l <- read.table('lg21_com.ld.gz',header=T)
l <- mutate(l,w1=floor(POS1/50000))
l <- mutate(l,w2=floor(POS2/50000))
l <- mutate(l,p=paste(w1,w2,sep="_"))
d <- l %>% group_by(p) %>% summarise(percentile=quantile(R.2,0.98,na.rm=T))
d1<-separate(d,p,into=c("w1","w2"),sep="_")
d1$w1 <- as.numeric(d1$w1)
d1$w2 <- as.numeric(d1$w2)
d1 <- arrange(d1,w1,w2)
write.table(d1,file="lg21w_com",quote=F,row.names=F,sep="\t")
##plot window LD results
#d1$w1 <- d1$w1*50000/1000000
#d1$w2 <- d1$w2*50000/1000000
#p<-ggplot(d,aes(x=w1,y=w2)) +
#     theme_bw() +
#     geom_tile(aes(fill=percentile)) + scale_fill_continuous(type = "viridis" )+
#     labs(title="LG26",x="Chromosome 26 (Mbp)",y="Chromosome 26 (Mbp)") +
#     theme(plot.title = element_text(size=30,hjust = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
#ggsave(filename = "lg26.svg",plot = p)
