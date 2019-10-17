library('data.table')
library('tidyverse')
library('stringr')
df<-fread('/path/to/fulldataframe')
summarydf<-df%>%group_by(chr, pos, pos2, ref, alt, gene_name, frameshift,stopgain)%>%summarise(zeros=n_distinct(sample[genotype=='0|0']),
                                                                                               het=n_distinct(sample[genotype=='1|0' | genotype=='0|1']),
                                                                                               ones=n_distinct(sample[genotype=='1|1']))
summarydf<-summarydf%>%mutate(raf=(het+2*min(zeros, ones))/4402)
overalldf<-as.data.frame(matrix(ncol=45, nrow=0))
                                        #overalldf[,c(4:6, 13)]<-as.character(overalldf[,c(4:6, 13)])
names(overalldf)<-c("chr", "pos", "pos2", "ref", "alt", "gene_name", "frameshift", "stopgain", "zeros", "het", "ones", "raf", "id", "gnomad_all_genome_af", "gnomad_afr_genome_af", "gnomad_amr_genome_af", "gnomad_asj_genome_af", "gnomad_eas_genome_af", "gnomad_fin_genome_af", "gnomad_nfe_genome_af", "gnomad_oth_genome_af", "gnomad_all_genome_ac", "gnomad_afr_genome_ac", "gnomad_amr_genome_ac", "gnomad_asj_genome_ac", "gnomad_eas_genome_ac", "gnomad_fin_genome_ac", "gnomad_nfe_genome_ac", "gnomad_oth_genome_ac", "gnomad_all_exome_af", "gnomad_afr_exome_af", "gnomad_amr_exome_af", "gnomad_asj_exome_af", "gnomad_eas_exome_af", "gnomad_fin_exome_af", "gnomad_nfe_exome_af", "gnomad_oth_exome_af", "gnomad_all_exome_ac", "gnomad_afr_exome_ac", "gnomad_amr_exome_ac", "gnomad_asj_exome_ac", "gnomad_eas_exome_ac", "gnomad_fin_exome_ac", "gnomad_nfe_exome_ac", "gnomad_oth_exome_ac")
names(overalldf)<-c('chr' ,'pos', 'pos2', 'ref', 'alt', 'gene_name', 'frameshift', 'stopgain', 'zeros', 'het', 'ones', 'raf', 'id', 'gnomad_all_genome_af', 'gnomad_nfe_genome_af', 'gnomad_all_exome_af', 'gnomad_nfe_exome_af')
for (i in 1:22){
name1=paste0('path/to/gnomad/nocomments/chr', i, '.PASS.post.beagle.4.1.gt.final.annot.first8.nocomments.vcf')
name2=paste0('/path/to/nocomments/chr', i, '.PASS.post.beagle.4.1.gt.final.exome.annot.first8.nocomments.vcf')
genomedf<-fread(name1)
exomedf<-fread(name2)

#GENOME AFs
genomedf$gnomad_all_genome_af<-as.numeric(str_match(genomedf$V8, "GNOMAD_ALL_AF=(.*?);")[,2])
genomedf$gnomad_afr_genome_af<-as.numeric(str_match(genomedf$V8, "GNOMAD_AFR_AF=(.*?);")[,2])
genomedf$gnomad_amr_genome_af<-as.numeric(str_match(genomedf$V8, "GNOMAD_AMR_AF=(.*?);")[,2])
genomedf$gnomad_asj_genome_af<-as.numeric(str_match(genomedf$V8, "GNOMAD_ASJ_AF=(.*?);")[,2])
genomedf$gnomad_eas_genome_af<-as.numeric(str_match(genomedf$V8, "GNOMAD_EAS_AF=(.*?);")[,2])
genomedf$gnomad_fin_genome_af<-as.numeric(str_match(genomedf$V8, "GNOMAD_FIN_AF=(.*?);")[,2])
genomedf$gnomad_nfe_genome_af<-as.numeric(str_match(genomedf$V8, "GNOMAD_NFE_AF=(.*?);")[,2])
genomedf$gnomad_oth_genome_af<-as.numeric(str_match(paste0(genomedf$V8, ";"), "GNOMAD_OTH_AF=(.*?);")[,2])

#GENOME ACs
genomedf$gnomad_all_genome_ac<-as.numeric(str_match(genomedf$V8, "GNOMAD_ALL_AC=(.*?);")[,2])
genomedf$gnomad_afr_genome_ac<-as.numeric(str_match(genomedf$V8, "GNOMAD_AFR_AC=(.*?);")[,2])
genomedf$gnomad_amr_genome_ac<-as.numeric(str_match(genomedf$V8, "GNOMAD_AMR_AC=(.*?);")[,2])
genomedf$gnomad_asj_genome_ac<-as.numeric(str_match(genomedf$V8, "GNOMAD_ASJ_AC=(.*?);")[,2])
genomedf$gnomad_eas_genome_ac<-as.numeric(str_match(genomedf$V8, "GNOMAD_EAS_AC=(.*?);")[,2])
genomedf$gnomad_fin_genome_ac<-as.numeric(str_match(genomedf$V8, "GNOMAD_FIN_AC=(.*?);")[,2])
genomedf$gnomad_nfe_genome_ac<-as.numeric(str_match(genomedf$V8, "GNOMAD_NFE_AC=(.*?);")[,2])
genomedf$gnomad_oth_genome_ac<-as.numeric(str_match(genomedf$V8, "GNOMAD_OTH_AC=(.*?);")[,2])
gnodf<-genomedf[,c(1:5,9:24)]
names(gnodf)[1:5]<-c('chr','pos','id','ref','alt')


#EXOME AFs
exomedf$gnomad_all_exome_af<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_ALL_AF=(.*?);")[,2])
exomedf$gnomad_all_exome_af<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_ALL_AF=(.*?);")[,2])
exomedf$gnomad_afr_exome_af<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_AFR_AF=(.*?);")[,2])
exomedf$gnomad_amr_exome_af<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_AMR_AF=(.*?);")[,2])
exomedf$gnomad_asj_exome_af<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_ASJ_AF=(.*?);")[,2])
exomedf$gnomad_eas_exome_af<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_EAS_AF=(.*?);")[,2])
exomedf$gnomad_fin_exome_af<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_FIN_AF=(.*?);")[,2])
exomedf$gnomad_nfe_exome_af<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_NFE_AF=(.*?);")[,2])
exomedf$gnomad_oth_exome_af<-as.numeric(str_match(paste0(exomedf$V8, ";"), "GNOMAD_EXOME_OTH_AF=(.*?);")[,2])

#EXOME ACs
exomedf$gnomad_all_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_ALL_AC=(.*?);")[,2])
exomedf$gnomad_all_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_ALL_AC=(.*?);")[,2])
exomedf$gnomad_afr_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_AFR_AC=(.*?);")[,2])
exomedf$gnomad_amr_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_AMR_AC=(.*?);")[,2])
exomedf$gnomad_asj_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_ASJ_AC=(.*?);")[,2])
exomedf$gnomad_eas_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_EAS_AC=(.*?);")[,2])
exomedf$gnomad_fin_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_FIN_AC=(.*?);")[,2])
exomedf$gnomad_nfe_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_NFE_AC=(.*?);")[,2])
exomedf$gnomad_oth_exome_ac<-as.numeric(str_match(exomedf$V8, "GNOMAD_EXOME_OTH_AC=(.*?);")[,2])

gnoexdf<-exomedf[,c(1:5,9:24)]
names(gnoexdf)[1:5]<-c('chr','pos','id','ref','alt')

joined_df1<-left_join(summarydf%>%filter(chr==i), gnodf)
joined_df2<-left_join(summarydf%>%filter(chr==i), gnoexdf)
joined_df<-left_join(joined_df1,joined_df2)
write.table(joined_df, file=paste0("joined_data_", i, '.test'))
overalldf<-rbind.data.frame(overalldf,joined_df)
}
write.table(overalldf, sep="\t", file = "joined_raf_dataframe", append = FALSE, quote = FALSE,row.names = FALSE, col.names=TRUE)
