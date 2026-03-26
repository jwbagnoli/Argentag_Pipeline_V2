#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

bam_path<-args[1]
yaml_path<-args[2]
Csel_type<- args[3]

##### Attention!!!!!!!!!!!!!!!!!! Install bambu via: https://github.com/GoekeLab/bambu/tree/Multiplex_Major_Patch , other version do not work!!!
library(bambu, quietly = T)
library(NanoporeRNASeq, quietly = T)
library(yaml, quietly = T)

ym<-read_yaml(yaml_path)

sample<- ym$sampleinfo$sample
outdir<- ym$general$outdir
ncores<- ym$general$nthreads
gtf_path<-ym$reference$gtf_path
fa_path<-ym$reference$genome_path

annotation <- prepareAnnotations(gtf_path)

dir.create(paste0(outdir, "/Bambu_",Csel_type ))

print(".....Creating read class file.")
readClassFile <- bambu(reads = bam_path, annotations = annotation, genome = fa_path, 
                       ncore = ncores, discovery = FALSE, quant = FALSE, 
                       demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE, 
                       yieldSize = 10000000, 
                       sampleNames = sample, cleanReads = TRUE,dedupUMI = TRUE, stranded = TRUE, trackReads = T)

saveRDS(readClassFile, paste0(outdir, "/Bambu_", Csel_type, "/", sample, "_readClassFile_", Csel_type, ".rds"))

print(".....Creating extended annotation.")
extendedAnno <- bambu(reads = readClassFile, annotations = annotation, genome = fa_path, 
                      ncore = ncores, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, 
                      verbose = FALSE, assignDist = FALSE, stranded = TRUE, trackReads = T)
saveRDS(extendedAnno, paste0(outdir,  "/Bambu_",Csel_type, "/", sample, "_extendedAnno_", Csel_type, ".rds"))

print("....Creating quant data")
quantData <- bambu(reads = readClassFile, annotations = extendedAnno, genome = fa_path, ncore = ncores, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, 
                   opt.em = list(degradationBias = FALSE), assignDist = TRUE, stranded = TRUE, trackReads = T)
saveRDS(quantData, paste0(outdir, "/Bambu_",Csel_type, "/", sample,  "_quantData_",Csel_type,".rds"))

print(".....Creating count matrices")
quantData.gene <- transcriptToGeneExpression(quantData[[1]])
counts <- assays(quantData.gene)$counts #selecting first sample

saveRDS(counts, paste0(outdir, "/Bambu_",Csel_type, "/", sample, "_counts_", Csel_type,".rds"))

rm(readClassFile,extendedAnno,quantData,quantData.gene, counts)
gc()