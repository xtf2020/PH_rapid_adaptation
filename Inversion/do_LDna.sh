#!/bin/bash
vcf="qc10k.vcf"
maf=0.1
for pop in all;
do
    out=$maf
    [ -f "$out/LDna.done" ] && continue # check
    rm -r $out
    :>$pop.LDna.cmd
    mkdir $out
    cut -f1 $pop >$out/sample.id
    for i in `seq 1 28`;
    do
        echo "vcftools --maf $maf --vcf $vcf --keep $out/sample.id --chr LG$i --geno-r2 -c | bgzip -@2 >$out/LG$i.ld.gz; \
		Rscript LDna.R $out/LG$i.ld.gz $out/LG$i.LDna $vcf $out/sample.id" >>$pop.LDna.cmd
    done
    ParaFly -c $pop.LDna.cmd -CPU 28
    # get PCA
    Rscript pca_and_boxplot_SeqVarTools.R $out
    touch $out/LDna.done
done

