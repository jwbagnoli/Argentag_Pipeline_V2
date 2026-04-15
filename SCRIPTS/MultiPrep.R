#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


yaml_path<- args[1]

### Libs
suppressPackageStartupMessages({library("ShortRead", quietly = T)
  library(dplyr, quietly = T)
  library(tidyr, quietly = T)
  library(pathviewr, quietly = T)
  library(ggplot2, quietly = T)
  library(ggridges, quietly = T)
  library(cowplot, quietly = T)
  library(yaml, quietly = T)
  library(foreach, quietly = T)
  library(doParallel, quietly = T)})

ym<-read_yaml(yaml_path)


nsamples<-length(ym$sampleinfo)
yams<-c()
snames<-c()
for (i in 1:nsamples){
  sname<-names(ym$sampleinfo)[i]
  if (ym$general$Start_stage == "taggy"){
    dir.create(paste0(ym$general$outdir, "/", sname))
  }
  sampleyam<-list()
  sampleyam$sampleinfo<-list()
  sampleyam$sampleinfo$sample_name<-sname
  sampleyam$sampleinfo$input_path<-as.character(ym$sampleinfo[i])
  sampleyam$general<-ym$general
  sampleyam$general$outdir<-paste0(ym$general$outdir, "/", sname)
  sampleyam$taggy<-ym$taggy
  sampleyam$filtering<-ym$filtering[[i]]
  sampleyam$reference<-ym$reference[[i]]
  sampleyam$tools<-ym$tools
  yaml::write_yaml(sampleyam, file = paste0(ym$general$outdir, "/", sname, "/", sname,".yaml"))
  yams<-c(yams,paste0(ym$general$outdir, "/", sname, "/", sname,".yaml") )
  snames<-c(snames, sname)
}
yams<-data.frame(samples = snames, yaml_path=yams)
write.table(yams,paste0(ym$general$outdir, "/yams.txt" ), row.names = F, col.names = F, quote = F, sep="\t")
