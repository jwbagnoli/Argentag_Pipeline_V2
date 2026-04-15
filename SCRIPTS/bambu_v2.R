#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

bam_path<-args[1]
yaml_path<-args[2]
Csel_name<- args[3]

##### Attention!!!!!!!!!!!!!!!!!! Install bambu via: https://github.com/GoekeLab/bambu/tree/Multiplex_Major_Patch , other version do not work!!!
suppressPackageStartupMessages({library(bambu, quietly = T)
library(NanoporeRNASeq, quietly = T)
library(yaml, quietly = T)
library(GGally, quietly = T)
library(tidyverse, quietly = T)
library(viridis)})
  
ym<-read_yaml(yaml_path)

sample<- ym$sampleinfo$sample_name
outdir<- ym$general$outdir
ncores<- ym$general$nthreads
gtf_path<-ym$reference$gtf_path
fa_path<-ym$reference$genome_path
Csel_type<-ym$filtering$Cellselection
gene_keys<- stringr::str_trim(c(unlist(strsplit(ym$reference$gene_keys, split = ",")), "BambuGene"))
species_names<- stringr::str_trim(c(unlist(strsplit(ym$reference$species_names, split = ",")), "Novel_Annotation"))

annotation <- prepareAnnotations(gtf_path)

dir.create(paste0(outdir, "/Bambu_",Csel_name ))

print(".....Creating read class file.")
readClassFile <- bambu(reads = bam_path, annotations = annotation, genome = fa_path, 
                       ncore = ncores, discovery = FALSE, quant = FALSE, 
                       demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE, 
                       yieldSize = 10000000, 
                       sampleNames = sample, cleanReads = TRUE,dedupUMI = TRUE, stranded = TRUE, trackReads = T)

saveRDS(readClassFile, paste0(outdir, "/Bambu_", Csel_name, "/", sample, "_readClassFile_", Csel_name, ".rds"))

print(".....Creating extended annotation.")
extendedAnno <- bambu(reads = readClassFile, annotations = annotation, genome = fa_path, 
                      ncore = ncores, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, 
                      verbose = FALSE, assignDist = FALSE, stranded = TRUE, trackReads = T)
saveRDS(extendedAnno, paste0(outdir,  "/Bambu_",Csel_name, "/", sample, "_extendedAnno_", Csel_name, ".rds"))

print("....Creating quant data")
quantData <- bambu(reads = readClassFile, annotations = extendedAnno, genome = fa_path, ncore = ncores, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, 
                   opt.em = list(degradationBias = FALSE), assignDist = TRUE, stranded = TRUE, trackReads = T)
saveRDS(quantData, paste0(outdir, "/Bambu_",Csel_name, "/", sample,  "_quantData_",Csel_name,".rds"))

print(".....Creating count matrices")
quantData.gene <- transcriptToGeneExpression(quantData[[1]])
counts <- assays(quantData.gene)$counts #selecting first sample

saveRDS(counts, paste0(outdir, "/Bambu_",Csel_name, "/", sample, "_counts_", Csel_name,".rds"))

# Create gene stats
#### split count matrix from bambu in different ones
split_counts_genes<- function(counts, patterns, names){
  counts_list<-list()
  genes<-rownames(counts)
  for (i in 1:length(patterns)){
    genes_tmp<-grep(x = genes, pattern = patterns[i], value = T)
    indcs = which(rownames(counts) %in% genes_tmp)
    counts_list[[i]] <- counts[indcs,]
  }
  names(counts_list)<-names
  return(counts_list)
}

counts_list<-split_counts_genes(counts = counts, patterns =gene_keys, names = species_names )

#### Create QC
qc_list<-list()
for (i in 1: length(counts_list)){
  qc_list[[i]]<-data.frame(BC =colnames(counts_list[[i]]),
                           UMIs= colSums(counts_list[[i]]),
                           Genes= colSums(counts_list[[i]]>0))
  
}  

names(qc_list)<-names(counts_list)
qc_df<-bind_rows(qc_list, .id = "Species")
qc_df$Species<- factor(qc_df$Species, levels = species_names)

#### Violin Plots
a<-ggplot()+
  geom_violin(data=qc_df, aes(x= Species, y=UMIs, fill=Species))+
  geom_jitter(data=qc_df, aes(x= Species, y=UMIs),alpha=0.2, size=1, width = 0.2)+
  theme_bw()+
  facet_wrap(~Species, scales="free", ncol=length(counts_list))+
  scale_fill_viridis(discrete=TRUE)+
  xlab("")+
  ggtitle("nUMIs")+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())+
  guides(fill="none")

b<-ggplot()+
  geom_violin(data=qc_df, aes(x= Species, y=Genes, fill=Species))+
  geom_jitter(data=qc_df, aes(x= Species, y=Genes),alpha=0.2, size=1, width = 0.2)+
  theme_bw()+
  facet_wrap(~Species, scales="free", ncol=length(counts_list))+
  scale_fill_viridis(discrete=TRUE)+
  xlab("")+
  ggtitle("nGenes")+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "bottom")+
  guides(fill="none")

qc<-cowplot::plot_grid(a,b, ncol=1, rel_heights = c(1,1))

ggsave(plot=qc, 
       filename = paste0(outdir,  "/Bambu_",Csel_name, "/", sample,"_Genes_UMIs_", Csel_name,".pdf" ),
       device = "pdf", 
       width = length(counts_list)*50, 
       height = 152,
       units= "mm"
)
ggsave(plot=qc, 
       filename = paste0(outdir,"/Bambu_",Csel_name, "/", sample,"_Genes_UMIs_", Csel_name,".jpg" ),
       device = "jpeg", 
       width = length(counts_list)*50, 
       height = 152,
       units= "mm"
)


### Correlation plots
qc_df_gene_wide<- qc_df[, c("BC","Species", "Genes")] %>% pivot_wider(names_from = Species, values_from = Genes) %>% column_to_rownames(var = "BC")
qc_df_UMI_wide<- qc_df[, c("BC","Species", "UMIs")] %>% pivot_wider(names_from = Species, values_from = UMIs) %>% column_to_rownames(var = "BC")
pm_genes <- ggpairs(qc_df_gene_wide)+theme_bw()+ggtitle("nGenes")
pm_UMIs <- ggpairs(qc_df_UMI_wide)+theme_bw()+ggtitle("nUMIs")


ggsave(plot=pm_genes, 
       filename = paste0(outdir, "/Bambu_",Csel_name, "/", sample,"_Genes_Corr_", Csel_name,".pdf" ),
       device = "pdf", 
       width = length(counts_list)*50, 
       height = (length(counts_list)*50)+10,
       units= "mm"
)
ggsave(plot=pm_genes, 
       filename = paste0(outdir, "/Bambu_",Csel_name, "/", sample,"_Genes_Corr_", Csel_name,".jpg" ),
       device = "jpeg", 
       width = length(counts_list)*50, 
       height = (length(counts_list)*50)+10,
       units= "mm"
)
ggsave(plot=pm_UMIs, 
       filename = paste0(outdir, "/Bambu_",Csel_name, "/", sample,"_UMIs_Corr_", Csel_name,".pdf" ),
       device = "pdf", 
       width = length(counts_list)*50, 
       height = (length(counts_list)*50)+10,
       units= "mm"
)
ggsave(plot=pm_UMIs, 
       filename = paste0(outdir, "/Bambu_",Csel_name, "/", sample,"_UMIs_Corr_", Csel_name,".jpg" ),
       device = "jpeg", 
       width = length(counts_list)*50, 
       height = (length(counts_list)*50)+10,
       units= "mm"
)

#### Clean
rm(readClassFile,extendedAnno,quantData,quantData.gene, counts)
gc()