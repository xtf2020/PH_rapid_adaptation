library(ggpubr)
library(magrittr)
library(dplyr)
library(ggforce)
library(magrittr)
library(gtools)

# read FWA.
fwa <- read.table('fresh.freq', header=T, comment.char="+")
fwa <- select(fwa,-HL)
names(fwa)[1] <- 'CHROM'
dat <- mutate(fwa, 
	 
	CHC=abs(HK-CH), 
	HZC=abs(HK-HZ),
	THC=abs(HK-TH), 
	Average=abs(HK-(TH+CH+HZ)/3)) %>%
	select(-HK,-TH,-CH,-HZ)
names(dat)[6:9] <- c('CH', 'HZ', 'TH', 'Average')
pop <- c(
"Chaohu Lake vs. Anadromous",
"Hongze Lake vs. Anadromous",
"Taihu Lake vs. Anadromous",
"Average")
names(pop) <-c ( 'CH', 'HZ','TH', 'Average')
abbr <- list(
CH='Chaohu Lake',
HZ='Hongze Lake',
TH='Taihu Lake',
HK='Anadromous'
)

hist_plt<-function(dat="", i="",j="", n=1) {

    p<-ggplot(dat, aes_string(i))+
    geom_histogram(breaks=seq(0,1,0.05))+
    scale_x_continuous( breaks = seq(0,1,0.1), limits=c(0,1))+
    labs(title=j, x="")+
    theme(
	text=element_text(family="Helvetica"),
	plot.title = element_text(hjust = 0.5,size=rel(1)),
	axis.text.x=element_text(color=rep(c("black","transparent"), 5)),
	panel.background = element_blank(),
	panel.border = element_blank(),
	axis.line = element_line(colour = "black"),
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	strip.background = element_blank()
	)
	if (n>1) p <- p + labs(y="")
    return(p)
}



plots_frq <- list()
plots_chg <- list()



# FWA
for (i in 1:length(abbr)) {
	ti <- abbr[[i]]
	po <- names(abbr)[i]
	plots_frq[[i]] <- hist_plt(fwa, po, ti, i)
}

# FWA change
for (i in 1:length(pop)) {
	ti <- pop[i]
	po <- names(pop)[i]
	plots_chg[[i]] <- hist_plt(dat, po, ti, i)
}

jpeg('FigS_fresh_change.jpg', height=500,width=1400)


fig1 <- ggarrange(plotlist=plots_frq, nrow=1)
fig1 <- annotate_figure(fig1, fig.lab = "A", fig.lab.size = 14, fig.lab.face = "bold")
fig2 <- ggarrange(plotlist=plots_chg, nrow=1)
fig2 <- annotate_figure(fig2, fig.lab = "B", fig.lab.size = 14, fig.lab.face = "bold")

#final fig.
ggarrange(fig1, fig2, ncol=1, title="ffe",heights=c(1,1))
dev.off()

