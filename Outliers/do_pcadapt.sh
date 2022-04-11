#!/bin/bash

for d in th/ hz/ ch/; 
do 
    cd $d; Rscript ../pcadapt.r; 
    cd ../; 
done

cat th/pcadapt.outliers.id hz/pcadapt.outliers.id \
ch/pcadapt.outliers.id \
| sort -V | uniq -c >outliers_pca.id

awk '$1==3 {print $2}' outliers_pca.id >fresh.pcadapt.outliers


