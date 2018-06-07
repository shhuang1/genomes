# create TAIR10 gene annotation files for TAIR10
library(biomaRt)
library(doMC)
library(GenomicFeatures)
library(plyr)
library(rtracklayer)

chr_alias = list('1'='chr1','2'='chr2','3'='chr3','4'='chr4','5'='chr5',
    'Mt'='chrM','Pt'='chrC')

write_feature_bed<-function(feature_gr,feature_name,output_file) {
    feature_df = as.data.frame(feature_gr,optional=T)
    if ('exon_rank' %in% colnames(feature_df)) {
        feature_rank = feature_df[,'exon_rank']
        name = feature_df[,'exon_name']
    } else {
        feature_rank = unlist(dlply(feature_df,'group_name',function(df) 1:nrow(df)))
        name = with(feature_df,paste0(group_name,'.',feature_name,'.',feature_rank))
    }
#    element_split = do.call("rbind",strsplit(feature_df[,'element'],"\\."))
#    element_name = element_split[,1]
#    element_rank = element_split[,2]
    feature_bed = with(feature_df,{
        data.frame(chrom = seqnames,
        start = start-1,
        end = end,
        name = paste0(group_name,"_",feature_name,"_",feature_rank),
        score = 0,
        strand = strand)
    })
    write.table(feature_bed,output_file,sep='\t',row.names=F,col.names=F,quote=F)
}

write_feature_bed_by_tx<-function(feature_gr,feature_name,output_file) {
    feature_df = as.data.frame(feature_gr,optional=T)
    if ('exon_rank' %in% colnames(feature_df)) {
        feature_df[,'feature_rank'] = feature_df[,'exon_rank']
    } else {
        feature_df = ddply(feature_df,'group_name',transform,feature_rank=order(seqnames,start,end))
    }
    feature_bed = with(feature_df,{
        data.frame(chrom = seqnames,
        start = start-1,
        end = end,
        name = paste0(group_name,"_",feature_name,"_",feature_rank),
        score = 0,
        strand = strand)
    })
    write.table(feature_bed,output_file,sep='\t',row.names=F,col.names=F,quote=F)
}

#txdb = makeTranscriptDbFromBiomart(biomart="plants_mart_22",dataset="athaliana_eg_gene",circ_seqs=c("Mt","Pt"))
#exons = exonsBy(txdb,"tx")
#write_feature_bed(exons,"exon",file.path(PROJ_DATA_PATH,"EnsemblPlants","release22","tair10.exons"))

write_feature_bed_all<-function(Gtf_dir) {
                               
    txdbGtf = makeTxDbFromGFF(file.path(Gtf_dir,"genes.gtf"),format="gtf")

                                        #rep_gene_models = read.table(file.path(PROJ_DATA_PATH,"TAIR10","Genes","TAIR10_genome_release","TAIR10_gene_lists","TAIR10_representative_gene_models"),stringsAsFactors=F,blank.lines.skip=T,skip=5)
                                        #txByGene = transcriptsBy(txdbGtf,"gene")
                                        #geneID = rep(names(txByGene),elementLengths(txByGene))
                                        #txByGeneDf = data.frame(geneID=geneID,txName=values(unlist(txByGene))[['tx_name']])
    
    exonsGtf = exonsBy(txdbGtf,"tx",use.names=T)
    write_feature_bed_by_tx(exonsGtf,"exon",file.path(Gtf_dir,'exons.bed'))
#write_feature_bed_by_tx(exonsGtf[rep_gene_models[,'V1']],"exon",file.path(PROJ_DATA_PATH,"TAIR10","Genes","TAIR10_genome_release","TAIR10_gff3","TAIR10_GFF3_genes_transposons.fixed2.chr.rep_genes.exons_by_tx.bed"))

    intronsGtf = intronsByTranscript(txdbGtf,use.names=T)
    write_feature_bed_by_tx(intronsGtf,"intron",file.path(Gtf_dir,'introns.bed'))
#write_feature_bed_by_tx(intronsGtf[rep_gene_models[,'V1']],"intron",file.path(PROJ_DATA_PATH,"TAIR10","Genes","TAIR10_genome_release","TAIR10_gff3","TAIR10_GFF3_genes_transposons.fixed2.chr.rep_genes.introns_by_tx.bed"))

    fiveUTRGtf = fiveUTRsByTranscript(txdbGtf,use.names=T)
    write_feature_bed_by_tx(fiveUTRGtf,"utr5",file.path(Gtf_dir,'utr5.bed'))
#write_feature_bed_by_tx(fiveUTRGtf[names(fiveUTRGtf) %in% rep_gene_models[,'V1']],"utr5",file.path(PROJ_DATA_PATH,"TAIR10","Genes","TAIR10_genome_release","TAIR10_gff3","TAIR10_GFF3_genes_transposons.fixed2.chr.rep_genes.utr5_by_tx.bed"))

    threeUTRGtf = threeUTRsByTranscript(txdbGtf,use.names=T)
    write_feature_bed_by_tx(threeUTRGtf,"utr3",file.path(Gtf_dir,'utr3.bed'))
#write_feature_bed_by_tx(threeUTRGtf[names(threeUTRGtf) %in% rep_gene_models[,'V1']],"utr3",file.path(PROJ_DATA_PATH,"TAIR10","Genes","TAIR10_genome_release","TAIR10_gff3","TAIR10_GFF3_genes_transposons.fixed2.chr.rep_genes.utr3_by_tx.bed"))

    genicGtf = genes(txdbGtf,single.strand.genes.only=FALSE)
    intergenic = gaps(unlist(genicGtf))
    
    genic2 = reduce(genicGtf,ignore.strand=TRUE)
    intergenic2 = gaps(unlist(genic2))
    intergenic2 = intergenic2[strand(intergenic2)=='*',]
    
    export(intergenic2,file.path(Gtf_dir,'intergenic.bed'),'BED')
}

write_feature_bed_all(Gtf_dir=file.path(PROJ_DATA_PATH,'Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2015-07-17-14-28-46/Genes'))

write_feature_bed_all(Gtf_dir=file.path(PROJ_DATA_PATH,'Illumina_iGenomesHomo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-2015-07-17-14-31-42/Genes/'))
