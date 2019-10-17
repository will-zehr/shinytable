library('data.table')
library('tidyverse')
##Needs vcfs to be in one directory, formatted as chr*.vcf, as 1:22. If you have a chr x, you should make an array of 1:22, X and use that for the for loop
df151<-read.csv('/path/to/genesofinterest', header=F)
names(df151)[1]<-'Gene_151'
pcdf<-fread('path/to/pcdf')
chrsummarydf<-as.data.frame(matrix(ncol=10,nrow=0))
names(chrsummarydf)<-c('CHR','LoF', 'NMD', 'FrameShift', 'Stop_Gain', 'Splice_Donor', 'Splice_Acceptor', 'Splice_Donor_or_Acceptor','NMD_Splice', 'Stop_Lost')
classsummarydf<-as.data.frame(matrix(0, ncol=10,nrow=3))
names(classsummarydf)<-c('Class','LoF', 'NMD', 'FrameShift', 'Stop_Gain', 'Splice_Donor', 'Splice_Acceptor', 'Splice_Donor_or_Acceptor','NMD_Splice', 'Stop_Lost')
genesummarydf<-as.data.frame(matrix(ncol=11,nrow=0))
names(genesummarydf)<-c('Gene','Chr','LoF', 'NMD', 'FrameShift', 'Stop.Gain', 'Splice.Donor', 'Splice.Acceptor', 'Splice.Donor.or.Acceptor','NMD.Splice', 'Stop.Lost')

for (i in 1:22){
    tablename<-paste0('chr', i, '.vcf')
    df<-fread(tablename)
    names(df)<-as.character(df[1,])
    df<-df[-1,]
    Extra<-df$Extra
    df<-df %>% separate(Extra, into=paste0("Extra", 1:2), "VARIANT_CLASS=")
    df<-df %>% separate(Extra2, into=c('Class', 'Extra3'), ';')
    df$Class<-as.factor(as.character(df$Class))
    df$Extra<-Extra
    df<-df%>%mutate(LoF.ind=grepl('HC', Extra),
                    NMD.ind=grepl('NMD', Consequence),
                    frameshift.ind=grepl('frameshift_variant', Consequence),
                    NMD.ind=grepl('NMD', Consequence),
                    SG.ind=grepl('stop_gain', Consequence),
                    sd.ind=grepl('splice_donor', Consequence),
                    sa.ind=grepl('splice_acceptor', Consequence),
                    sl.ind=grepl('stop_lost', Consequence),
                    SE.flag=grepl('SINGLE_EXON', Extra),
                    NN.flag=grepl('NAGNAG', Extra),
                    PHYS.weak.flag=grepl('PHYLOCSF_WEAK', Extra),
                    PHYS.UO.flag=grepl('PHYLOCSF_UNLIKELY_ORF', Extra),
                    NCS.flag=grepl('NON_CAN_SPLICE', Extra),
                    ET.filter=grepl('END_TRUNC', Extra),
                    I.CDS.filter=grepl('INCOMPLETE_CDS', Extra),
                    Exon.intron.filter=grepl('EXON_INTRON_UNDEF', Extra),
                    Small.intron.filter=grepl('SMALL_INTRON', Extra),
                    Anc.filter=grepl('ANC_ALLELE', Extra),
                    NDD.filter=grepl('NON_DONOR_DISRUPTING', Extra),
                    NAD.filter=grepl('NON_ACCEPTOR_DISRUPTING', Extra),
                    Rescue.donor.filter=grepl('RESCUE_DONOR', Extra),
                    Rescue.acceptor.filter=grepl('RESCUE_ACCEPTOR', Extra),
                    GCGT.filter=grepl('GC_TO_GT_DONOR', Extra),
                    GCGT.filter=grepl('UTR_SPLICE', Extra)
                    )%>%mutate(
                            sd.sa=sa.ind+sd.ind,
                            sd.sa.NMD=sa.ind*NMD.ind+sd.ind*NMD.ind,
                            filter.ind=PHYS.weak.flag+PHYS.UO.flag+NCS.flag+
                                ET.filter+I.CDS.filter+Exon.intron.filter+Small.intron.filter+
                                Anc.filter+NDD.filter+NAD.filter+Rescue.donor.filter+Rescue.acceptor.filter+
                                GCGT.filter+GCGT.filter
                        )
    genedf<-left_join(df, pcdf, by=c('Gene'='V4'))
    genedf<-genedf%>%mutate(ind_151=ifelse(V5 %in% df151$Gene_151, 1,0)) %>% separate(Location, into=c('Chr', 'Pos'), ":") %>% separate(Pos, into=c('Start','End'), '-')
    genedf<-genedf%>%mutate(Start=as.numeric(Start),End=ifelse(is.na(End), as.numeric(Start)+1, End))
    genedf2<-cbind(genedf[,2:ncol(genedf)],genedf[1])
    tablename2<-paste0('chr_gene', i, '.withfilters.vcf')
    write.table(genedf2, sep="\t", file = tablename2, append = FALSE, quote = FALSE,row.names = FALSE, col.names=TRUE)}
