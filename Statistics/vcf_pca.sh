in=$1
out=$2
pop=$3
title=$4
pop=${pop:-2}
if [ $# -lt 2 ] || [ $1 = -h ];then
    echo "$0 [vcf] [out] [popmap] [title]"
    exit 0
fi
#plink

plink --vcf $in --out out --make-bed --allow-extra-chr --keep-allele-order --chr-set 28

if [ -f $pop ];then
    awk '{print $1"\t"$2}' out.fam | paste - $pop >fam.id
    
else
    awk -v p=$pop '{n=substr($2,1,p);print $1,$2,$1,n}' out.fam >fam.id
fi

plink --bfile out --out out --make-grm-bin --allow-extra-chr --chr-set 28 --update-ids fam.id
#gcta64
gcta64 --grm out --pca --out pca
rm out*
#ggplot2
R --vanilla <<eof
library('ggplot2')
dat<-read.table("pca.eigenvec",header=F)
var<-read.table("pca.eigenval",header=F)
names(dat)[2:4] <- c("Populations","PC1","PC2")
lab1<-paste("PC1"," (",sprintf("%.2f%%", 100*var[1,]/sum(var)),")", sep="")
lab2<-paste("PC2"," (",sprintf("%.2f%%", 100*var[2,]/sum(var)),")", sep="")
pdf("$out.pdf")
ggplot(dat,aes(x=PC1,y=PC2,group=Populations,colour=Populations,shape=Populations))+
labs(title="$title",x=lab1,y=lab2)+
geom_point()+
stat_ellipse()+
theme(plot.title=element_text(hjust=0.5))+
scale_shape_manual(values=seq(0,15))+
theme_classic()
dev.off()
eof
