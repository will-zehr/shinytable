Commands & R script to make a shiny table for collaborators to use in Norway. Check it out! 

[woverton.shinyapps.io/filtertable2/](https://woverton.shinyapps.io/filtertable2/)


# HUNT LoF Annotation

### Goal: annotate variants using whole genome sequence data from 2201 samples in HUNT

##  Analysis overview

1) annotate variants using loftee

```console
for f in `ls /path/to/vcfs`; do base=`basename $f .vcf`; echo "path/to/vep -i $f  --offline --plugin LoF,loftee_path:path/to.vep/Plugins/,human_ancestor_fa:/path/to/.vep/Plugins/human_ancestor.fa.gz,conservation_file:/path/to/.vep/Plugins/phylocsf_gerp.sql  -o $base.vcf --force_overwrite --cache --dir /path/to/.vep/ --dir_plugins /path/to/.vep/Plugins/"; done > submit_VEP.sh
```

2) parse variants using rscript filterparse.R, with output chr_gene*.withfilters.vcf

3) filter variants so we only wind up with quality HCLOF variants

```console
for i in `seq 1 22`; do awk '{ if ($19 == "TRUE" && $43==0) { print } }' chr_gene$i.withfilters.vcf > chr$i.filtered.vcf; done
mv *filtered.vcf FilteredVariants
```

4) connect sample genotype data to the variant data, using Brooke's parseBCFtoolsQuery.py script

```console
for i in `seq 1 22`; do for f in `ls path/to/vcfs`; do bcftools query -R chr$i.filtered.vcf -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' $f | python /path/to/parseBCFtoolsQuery.py -f - >> chr$i.filtered_sample_genotypes.dat; done; done
```

5) get summary data (counts of 1|1, 0|1, and 0|0) per snp, and compare to gnomad
```R
R
library('data.table')
library('tidyverse')
df2<-fread('all_chr_variants_filtered.vcf')
variants<-select(df2, Chr, Start, End, Allele, Gene, Feature, frameshift.ind, SG.ind, V5)%>%rename(gene_name=V5)
distvariants<-unique(variants)%>%group_by(Chr, Start, End, Allele, Gene, gene_name)%>%summarise(frameshift=as.factor(max(frameshift.ind)),stopgain=as.factor(max(SG.ind)))
write.table(distvariants, file='all_variants_unique_filtered.bed', sep="\t", append = FALSE, quote = FALSE,row.names = FALSE, col.names=FALSE)
```
```console
for i in `seq 1 22`; do awk 'NR > 0 {print $0}' chr$i.filtered_sample_genotypes.dat >> all_chr_genotypes.filtered.dat; done
```

6) run r script joinscript.R to output fulldataframe, which contains count information for all variants

7) get gnomad allele frequencies at our variant locations
```console
python gnomadAnnotate.py
```

8) run rscript gnomadjoinscript.R to join our AF info with gnomad AF info

And the app is ready to be made. See [code here](./app.R)
