#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


yaml_path<- args[1]


### Libs
library("ShortRead", quietly = T)
library(dplyr, quietly = T)
library(tidyr, quietly = T)
library(pathviewr, quietly = T)
library(ggplot2, quietly = T)
library(ggridges, quietly = T)
library(cowplot, quietly = T)
library(yaml, quietly = T)
library(foreach, quietly = T)
library(doParallel, quietly = T)

ym<-read_yaml(yaml_path)

outdir<- ym$general$outdir
chunks<- ym$general$nchunks
ncores<- ym$general$nthreads

automated_sel<-ym$filtering$automated_sel
ncells_sel<-ym$filtering$ncells_sel
custom_sel<-ym$filtering$custom_sel
ncells<-ym$filtering$ncells
thr.minreads<-ym$filtering$minreads
thr.trimlength.mean<-ym$filtering$minreadlength
registerDoParallel(cores=ncores)


# do in parallel cause vfc
fq_files<- list.files(paste0(outdir, "/taggy_demux/fastq"), full.names = T, pattern = "*[[:digit:]].fastq.gz")

chunk_list<-split(fq_files,ceiling(seq_along(fq_files) / (96/chunks)))

dir.create(paste0(outdir, "/Filtering/"))
dir.create(paste0(outdir, "/Filtering/Stats_tmp/"))
for (i in 1:length(chunk_list)) {
  dir.create(paste0(outdir,"/Filtering/Stats_tmp/Chunk_", i))
}

foreach(j=1:chunks) %dopar% { 
  sink(paste0(outdir, "/Filtering/Stats_tmp/StatsChecklog.txt"), append=TRUE) #open sink file and add output
  for (k in 1:length(chunk_list[[j]])){
    file_name<-as.character(strsplit(last(strsplit(chunk_list[[j]][k], split ="/", fixed=T)[[1]]), split =".", fixed=T)[[1]][1])
    cat(paste0("Processing file ", file_name , " in chunk ", j, ".\n"))
    fastQ <- readFastq(chunk_list[[j]][k])
    df<-as.data.frame(sread(fastQ) %>% as.character())
    colnames(df)<-"trimmed_seq"
    rownames(df)<-ShortRead::id(fastQ) %>% as.character()
    df$ReadID<- sapply(rownames(df), function(x){strsplit(x, split="#")[[1]][2]})
    df$BC_Triplet<- sapply(rownames(df), function(x){strsplit(x, split="_")[[1]][1]})
    df$UMI_seq<- sapply(rownames(df), function(x){strsplit(strsplit(x, split="_")[[1]][2], split="#")[[1]][1]})
    df$trimmed_length<- nchar(df$trimmed_seq)
    saveRDS(df, paste0(outdir, "/Filtering/Stats_tmp/Chunk_", j, "/File_",file_name,".rds"))
    # Create stats
    df_triplet<- df %>% group_by(BC_Triplet) %>% summarise(reads=n(), umis= length(unique(UMI_seq)), 
                                                           mean_length_trimmed=mean(trimmed_length),medium_length_trimmed=median(trimmed_length))
    
    rm(df)
    saveRDS(df_triplet, paste0(outdir, "/Filtering/Stats_tmp/Chunk_", j, "/File_",file_name,"_perBCstats.rds"))
    rm(df_triplet)
    gc()
  }
}

stats_list<-list()
for (i in 1:chunks){
  stats_files<-list.files(paste0(outdir, "/Filtering/Stats_tmp/Chunk_", i), full.names = T, pattern = "*._perBCstats.rds")
  stats_list[[i]]<-list()
  for (j in 1:length(stats_files)){
    stats_list[[i]][[j]]<-readRDS(stats_files[j])
  }
}


stats_list2<-list()
for (i in 1:chunks){
  stats_list2[[i]]<-bind_rows(stats_list[i])
}

stats<- bind_rows(stats_list2)
saveRDS(stats, paste0(outdir, "/Filtering/perBC_stats.rds"))

rm(stats_list,stats_list2)

# selection of cells 
### using automated cell selection
if (isTRUE(automated_sel)){
  readys_reads<-data.frame(reads=sort(stats$reads, decreasing = T))
  readys_reads$index<-as.numeric(rownames(readys_reads))
  thr.index_reads<-as.numeric(find_curve_elbow(readys_reads[1:50000,c("index", "reads")], plot_curve = F))
  thr.reads<-as.numeric(readys_reads$reads[thr.index_reads])
  
  if (isEmpty(thr.index_reads)){
    thr.index_reads<-50000
    print("Warning: Could not set elbow based on reads. taking 50000 Barcodes.")
  }
  
  a<-ggplot()+
    geom_point(data=readys_reads[1:50000,], aes(x=index, y=reads))+
    geom_text(aes(x=thr.index_reads, y=thr.reads, label=paste0(thr.index_reads, " cells")), hjust=-0.2)+
    geom_vline(xintercept = thr.index_reads)+
    scale_y_log10()+
    theme_bw()+
    ggtitle(paste0(thr.reads, " reads"))
  
  readys_length<-data.frame(length=sort(filter(stats, reads >= thr.reads)$mean_length_trimmed, decreasing = T))
  readys_length$index<-as.numeric(rownames(readys_length))
  
  thr.index_length<-as.numeric(find_curve_elbow(readys_length[c("index", "length")], plot_curve = T))
  thr.length<-floor(as.numeric(readys_length$length[thr.index_length]))
  
  if (isEmpty(thr.index_length)){
    thr.length<-mean_length_trimmed
    thr.index_length<-as.numeric(max(readys_length$index[thr.length]))
    print("Warning: Could not set elbow based on mean read length. taking provided manual read length.")
  }
  
  b<-ggplot()+
    geom_point(data=readys_length, aes(x=index, y=length))+
    geom_text(aes(x=thr.index_length, y=thr.length, label=paste0(thr.index_length, " bp")), hjust=-0.2)+
    geom_vline(xintercept = thr.index_length)+             
    theme_bw()+
    ggtitle(paste0(thr.length, " bp"))
  
  elbow<-cowplot::plot_grid(a,b, ncol = 2)
  elbow
  
  ggsave(plot=elbow, 
         filename = paste0(outdir, "/Filtering/automated_Cellselection.pdf" ),
         device = "pdf", 
         width = 310, 
         height = 152,
         units= "mm"
  )
  ggsave(plot=elbow, 
         filename = paste0(outdir, "/Filtering/automated_Cellselection.jpg" ),
         device = "jpeg", 
         width = 310, 
         height = 152,
         units= "mm"
  )
  
  stats_filt_auto<-filter(stats,reads >= thr.reads  & mean_length_trimmed > thr.length)
  write.table(stats_filt_auto$BC_Triplet,paste0(outdir, "/Filtering/Retained_BCs_automated.tsv"), quote = F, col.names = F, row.names = F)
  
  if (nrow(stats_filt_auto) == 0){
    print("Warning: Automated cell selection failed.")
  }
  
}

### using ncells cell selection
if (isTRUE(ncells_sel)){
  readys_ncells<-data.frame(reads=sort(filter(stats,mean_length_trimmed >= thr.trimlength.mean )$reads, decreasing = T))
  readys_ncells$index<-as.numeric(rownames(readys_ncells))
  if (ncells > nrow(readys_ncells)){
    ncells <- as.numeric(nrow(readys_ncells))
    print("Warning: number of targeted cells higher than number of Barcodes above read length threshold. Taking all above threshold for ncells BC selection")
  }
  write.table(ncells, paste0(outdir, "/Filtering/updated_ncells.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
  thr.reads_elbow_ncells<-as.numeric(readys_ncells$reads[ncells])
  
  elbow_ncells<-ggplot()+
    geom_point(data=readys_ncells[1:50000,], aes(x=index, y=reads), size=0.2, colour="grey")+
    geom_text(aes(x=ncells, y=thr.reads_elbow_ncells, label=paste0(ncells, " cells")), hjust=-0.2)+
    geom_vline(xintercept = ncells)+
    scale_y_log10()+
    theme_bw()+
    ggtitle(paste0(thr.reads_elbow_ncells, " reads"))
  elbow_ncells
  
  
  
  ggsave(plot=elbow_ncells, 
         filename = paste0(outdir, "/Filtering/",ncells ,"cells_Cellselection.pdf" ),
         device = "pdf", 
         width = 155, 
         height = 152,
         units= "mm"
  )
  ggsave(plot=elbow_ncells, 
         filename = paste0(outdir, "/Filtering/",ncells ,"cells_Cellselection.jpg" ),
         device = "jpeg", 
         width = 155, 
         height = 152,
         units= "mm"
  )
  
  rm(elbow, elbow_ncells)
  
  stats_filt_ncells<-filter(stats,reads >=  thr.reads_elbow_ncells & mean_length_trimmed > thr.trimlength.mean)
  write.table(stats_filt_ncells$BC_Triplet,paste0(outdir, "/Filtering/Retained_BCs_",ncells ,"cells.tsv"), quote = F, col.names = F, row.names = F)     
}


### using custom cell selection
if (isTRUE(custom_sel)){
  stats_filt_custom<-filter(stats,reads >=  thr.minreads & mean_length_trimmed > thr.trimlength.mean)
  
  if (nrow(stats_filt_custom)  == 0){
      print("Warning: No barcodes above minimum number of reads and readlength threshold. custom cell selection failed.")
  }

  
  
  write.table(stats_filt_custom$BC_Triplet,paste0(outdir, "/Filtering/Retained_BCs_custom.tsv"), quote = F, col.names = F, row.names = F)     
}

### global stats
nCselection<- sum(automated_sel, ncells_sel,custom_sel )
### density
plot_df<-list()
plot_df$all_BCs<-stats

if (isTRUE(automated_sel)){
  plot_df$automated<-stats_filt_auto
}
if (isTRUE(ncells_sel)){
  plot_df[[paste0(ncells,"cells")]]<-stats_filt_ncells
}
if (isTRUE(custom_sel)){
  plot_df$custom<-stats_filt_custom
}

plot_df<- bind_rows(plot_df, .id = "Cell_selection")
plot_df$Cell_selection<-factor(plot_df$Cell_selection, levels= c("all_BCs", "automated", paste0(ncells,"cells"), "custom"))

plot_df$Cell_selection<-droplevels(plot_df$Cell_selection)
# if (Csel_type =="both"){
#   plot_df<-list(stats,stats_filt,stats_filt_ncells  )
#   names(plot_df)<-c("all_BCs", "automated", paste0(ncells,"cells"))
#   plot_df<- bind_rows(plot_df, .id = "Cell_selection")
#   plot_df$Cell_selection<-factor(plot_df$Cell_selection, levels= c("all_BCs", "automated", paste0(ncells,"cells")))
# }else if(Csel_type =="auto"){
#   plot_df<-list(stats,stats_filt  )
#   names(plot_df)<-c("all_BCs", "automated")
#   plot_df<- bind_rows(plot_df, .id = "Cell_selection")
#   plot_df$Cell_selection<-factor(plot_df$Cell_selection, levels= c("all_BCs", "automated"))
# }else{
#   plot_df<-list(stats,stats_filt_ncells  )
#   names(plot_df)<-c("all_BCs", paste0(ncells,"cells"))
#   plot_df<- bind_rows(plot_df, .id = "Cell_selection")
#   plot_df$Cell_selection<-factor(plot_df$Cell_selection, levels= c("all_BCs",  paste0(ncells,"cells")))
# }

cols<- c("grey", "lightgreen", "lightblue", "purple4")
names(cols)<- c("all_BCs", "automated", paste0(ncells,"cells"), "custom")
A<-ggplot()+
  geom_density_ridges(data=plot_df, aes(x=reads,y=Cell_selection, fill=Cell_selection), alpha=0.5)+
  scale_x_log10()+
  theme_bw()+
  ylab("")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "none")
B<-ggplot()+
  geom_density_ridges(data=plot_df, aes(x=umis,y=Cell_selection, fill=Cell_selection), alpha=0.5)+
  scale_x_log10()+
  theme_bw()+
  ylab("")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "none")
C<-ggplot()+
  geom_density_ridges(data=plot_df, aes(x=mean_length_trimmed,y=Cell_selection, fill=Cell_selection), alpha=0.5)+
  scale_x_log10()+
  theme_bw()+
  ylab("")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "none")
D<-ggplot()+
  geom_density_ridges(data=plot_df, aes(x=medium_length_trimmed,y=Cell_selection, fill=Cell_selection), alpha=0.5)+
  scale_x_log10()+
  theme_bw()+
  ylab("")+
  scale_fill_manual(values =cols)+
  theme(legend.position = "none")

Reads_UMIs_Length_stats_density<-plot_grid(A,C,B,D,ncol = 2)
Reads_UMIs_Length_stats_density

ggsave(plot=Reads_UMIs_Length_stats_density, 
       filename = paste0(outdir, "/Filtering/stats_density.pdf" ),
       device = "pdf", 
       width = 200, 
       height = 50+(50*nCselection),
       units= "mm"
)
ggsave(plot=Reads_UMIs_Length_stats_density, 
       filename = paste0(outdir, "/Filtering/stats_density.jpeg" ),
       device ="jpeg", 
       width = 200, 
       height = 50+(50*nCselection),
       units= "mm"
)

rm(Reads_UMIs_Length_stats_density )



### dotplot
stats_plot<-sample_n(stats, size = 100000, replace = F)
plot_list<-list()
#### automated
if (isTRUE(automated_sel)){
  stats_left_up<-filter(stats_plot,reads <  thr.reads & mean_length_trimmed > thr.length)
  stats_left_down<-filter(stats_plot,reads <  thr.reads & mean_length_trimmed < thr.length)
  stats_right_up<-filter(stats_plot,reads >=  thr.reads & mean_length_trimmed > thr.length)
  stats_right_down<-filter(stats_plot,reads >=  thr.reads & mean_length_trimmed < thr.length)
  
  n_left_up<-as.numeric(nrow(filter(stats,reads <  thr.reads & mean_length_trimmed > thr.length)))
  n_left_down<-as.numeric(nrow(filter(stats,reads <  thr.reads & mean_length_trimmed < thr.length)))
  n_right_up<-as.numeric(nrow(filter(stats,reads >=  thr.reads & mean_length_trimmed > thr.length)))
  n_right_down<-as.numeric(nrow(filter(stats,reads >=  thr.reads & mean_length_trimmed < thr.length)))
  
  reads_trimlength_dot<-ggplot()+
    geom_point(data=stats_left_up, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_point(data=stats_left_down, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_point(data=stats_right_up, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="lightgreen")+
    geom_point(data=stats_right_down, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_text(aes(x=0.5*thr.reads, y=10000, label=n_left_up), colour="grey50", hjust = 0.5)+
    geom_text(aes(x=0.5*thr.reads, y=10, label=n_left_down), colour="grey50", hjust = 0.5)+
    geom_text(aes(x=thr.reads+(0.5*(max(stats_plot$reads)-thr.reads)), y=10000, label=n_right_up), colour="lightgreen", hjust = 0.5)+
    geom_text(aes(x=thr.reads+(0.5*(max(stats_plot$reads)-thr.reads)), y=10, label=n_right_down), colour="grey50", hjust = 0.5)+
    geom_vline(xintercept = thr.reads)+
    geom_hline(yintercept = thr.length)+
    scale_x_log10()+
    scale_y_log10()+
    theme_bw()+
    ggtitle("Automated Cell Selection", subtitle = paste0(">= ",thr.reads , " reads / >= ",thr.length, " mean trimmed length"))
  plot_list$automated<-reads_trimlength_dot
}


#### custom cell number
if (isTRUE(ncells_sel)){
  stats_left_up<-filter(stats_plot,reads <  thr.reads_elbow_ncells & mean_length_trimmed > thr.trimlength.mean)
  stats_left_down<-filter(stats_plot,reads <  thr.reads_elbow_ncells & mean_length_trimmed < thr.trimlength.mean)
  stats_right_up<-filter(stats_plot,reads >=  thr.reads_elbow_ncells & mean_length_trimmed > thr.trimlength.mean)
  stats_right_down<-filter(stats_plot,reads >=  thr.reads_elbow_ncells & mean_length_trimmed < thr.trimlength.mean)
  
  n_left_up_ncells<-as.numeric(nrow(filter(stats,reads <  thr.reads_elbow_ncells & mean_length_trimmed > thr.trimlength.mean)))
  n_left_down_ncells<-as.numeric(nrow(filter(stats,reads <  thr.reads_elbow_ncells & mean_length_trimmed < thr.trimlength.mean)))
  n_right_up_ncells<-as.numeric(nrow(filter(stats,reads >=  thr.reads_elbow_ncells & mean_length_trimmed > thr.trimlength.mean)))
  n_right_down_ncells<-as.numeric(nrow(filter(stats,reads >=  thr.reads_elbow_ncells & mean_length_trimmed < thr.trimlength.mean)))
  
  reads_trimlength_dot_ncells<-ggplot()+
    geom_point(data=stats_left_up, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_point(data=stats_left_down, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_point(data=stats_right_up, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="lightgreen")+
    geom_point(data=stats_right_down, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_text(aes(x=0.5*thr.reads_elbow_ncells, y=10000, label=n_left_up_ncells), colour="grey50", hjust = 0.5)+
    geom_text(aes(x=0.5*thr.reads_elbow_ncells, y=10, label=n_left_down_ncells), colour="grey50", hjust = 0.5)+
    geom_text(aes(x=thr.reads_elbow_ncells+(0.5*(max(stats_plot$reads)-thr.reads_elbow_ncells)), y=10000, label=n_right_up_ncells), colour="lightgreen", hjust = 0.5)+
    geom_text(aes(x=thr.reads_elbow_ncells+(0.5*(max(stats_plot$reads)-thr.reads_elbow_ncells)), y=10, label=n_right_down_ncells), colour="grey50", hjust = 0.5)+
    geom_vline(xintercept = thr.reads_elbow_ncells)+
    geom_hline(yintercept = thr.trimlength.mean)+
    scale_x_log10()+
    scale_y_log10()+
    theme_bw()+
    ggtitle(paste0(ncells, " cells Selection"), subtitle = paste0(">= ",thr.reads_elbow_ncells , " reads / >= ",thr.trimlength.mean, " mean trimmed length"))
  plot_list$ncells<-reads_trimlength_dot_ncells
}

#### custom reads
if (isTRUE(custom_sel)){
  stats_left_up<-filter(stats_plot,reads <  thr.minreads & mean_length_trimmed > thr.trimlength.mean)
  stats_left_down<-filter(stats_plot,reads <  thr.minreads & mean_length_trimmed < thr.trimlength.mean)
  stats_right_up<-filter(stats_plot,reads >=  thr.minreads & mean_length_trimmed > thr.trimlength.mean)
  stats_right_down<-filter(stats_plot,reads >=  thr.minreads & mean_length_trimmed < thr.trimlength.mean)
  
  n_left_up_cust<-as.numeric(nrow(filter(stats,reads <  thr.minreads & mean_length_trimmed > thr.trimlength.mean)))
  n_left_down_cust<-as.numeric(nrow(filter(stats,reads <  thr.minreads & mean_length_trimmed < thr.trimlength.mean)))
  n_right_up_cust<-as.numeric(nrow(filter(stats,reads >=  thr.minreads & mean_length_trimmed > thr.trimlength.mean)))
  n_right_down_cust<-as.numeric(nrow(filter(stats,reads >=  thr.minreads & mean_length_trimmed < thr.trimlength.mean)))
  
  reads_trimlength_dot_custom<-ggplot()+
    geom_point(data=stats_left_up, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_point(data=stats_left_down, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_point(data=stats_right_up, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="lightgreen")+
    geom_point(data=stats_right_down, aes(x=reads, y=mean_length_trimmed), size=0.5, alpha=0.5, colour="grey")+
    geom_text(aes(x=0.5*thr.minreads, y=10000, label=n_left_up_cust), colour="grey50", hjust = 0.5)+
    geom_text(aes(x=0.5*thr.minreads, y=10, label=n_left_down_cust), colour="grey50", hjust = 0.5)+
    geom_text(aes(x=thr.minreads+(0.5*(max(stats_plot$reads)-thr.minreads)), y=10000, label=n_right_up_cust), colour="lightgreen", hjust = 0.5)+
    geom_text(aes(x=thr.minreads+(0.5*(max(stats_plot$reads)-thr.minreads)), y=10, label=n_right_down_cust), colour="grey50", hjust = 0.5)+
    geom_vline(xintercept = thr.minreads)+
    geom_hline(yintercept = thr.trimlength.mean)+
    scale_x_log10()+
    scale_y_log10()+
    theme_bw()+
    ggtitle("Custom Cell Selection", subtitle = paste0(">= ",thr.minreads , " reads / >= ",thr.trimlength.mean, " mean trimmed length"))
  plot_list$custom<-reads_trimlength_dot_custom
}

width_plot<-150*nCselection

if (length(plot_list) == 1){
  Reads_Length_stats_dot<-plot_list[[1]]
}else{
  Reads_Length_stats_dot<-plot_grid(plotlist = plot_list, ncol=nCselection)
}

ggsave(plot=Reads_Length_stats_dot, 
       filename = paste0(outdir, "/Filtering/stats_dot.pdf" ),
       device = "pdf", 
       width = width_plot, 
       height = 152,
       units= "mm"
)
ggsave(plot=Reads_Length_stats_dot, 
       filename = paste0(outdir, "/Filtering/stats_dot.jpeg" ),
       device = "jpeg", 
       width = width_plot, 
       height = 152,
       units= "mm"
)
rm(stats_left_up,stats_left_down,stats_right_up,stats_right_down)




## select read IDs automated
chunk_list2<- list()
for (i in 1:chunks){
  chunk_list2[[i]]<- list.files(paste0(outdir, "/Filtering/Stats_tmp/Chunk_", i), full.names = T, pattern = "*[[:digit:]].rds")
}

foreach(j=1:chunks) %dopar% { 
  sink(paste0(outdir, "/Filtering/Stats_tmp/filterlog.txt"), append=TRUE) #open sink file and add output
  for (k in 1:length(chunk_list2[[j]])){
    file_name<-as.character(strsplit(strsplit(last(strsplit(chunk_list2[[j]][k], split ="/", fixed=T)[[1]]), split =".", fixed=T)[[1]][1], split="_")[[1]][2])
    cat(paste0("Processing file ", file_name , " in chunk ", j, ".\n"))
    # read in original fastq after grep subsetting
    df <- readRDS(chunk_list2[[j]][k])
    if (isTRUE(automated_sel)){
      df_filt <- filter(df, BC_Triplet %in% stats_filt_auto$BC_Triplet)
      if (nrow(df_filt) > 0){
        df_filt$ReadID_tagged <- paste0(df_filt$BC_Triplet, "_", df_filt$UMI_seq, "#", df_filt$ReadID)
        saveRDS(df_filt$ReadID_tagged, paste0(outdir,"/Filtering/Stats_tmp/Chunk_", j, "/File_",file_name,"_filtered_reads_automated.rds"))
      }
      rm(df_filt)
    }
    if (isTRUE(ncells_sel)){
      df_filt_ncells <- filter(df, BC_Triplet %in% stats_filt_ncells$BC_Triplet)
      if (nrow(df_filt_ncells) > 0){
        df_filt_ncells$ReadID_tagged <- paste0(df_filt_ncells$BC_Triplet, "_", df_filt_ncells$UMI_seq, "#", df_filt_ncells$ReadID)
        saveRDS(df_filt_ncells$ReadID_tagged, paste0(outdir,"/Filtering/Stats_tmp/Chunk_", j, "/File_",file_name,"_filtered_reads_",ncells ,"cells.rds"))
      }
      rm(df_filt_ncells)
    }
    if (isTRUE(custom_sel)){
      df_filt_custom <- filter(df, BC_Triplet %in% stats_filt_custom$BC_Triplet)
      if (nrow(df_filt_custom) > 0){
        df_filt_custom$ReadID_tagged <- paste0(df_filt_custom$BC_Triplet, "_", df_filt_custom$UMI_seq, "#", df_filt_custom$ReadID)
        saveRDS(df_filt_custom$ReadID_tagged, paste0(outdir,"/Filtering/Stats_tmp/Chunk_", j, "/File_",file_name,"_filtered_reads_custom.rds"))
      }
      rm(df_filt_custom)
    }
    rm(df)
  }
}


if (isTRUE(automated_sel)){
  filter_list<-list()
  for (i in 1:chunks){
    filter_files<-list.files(paste0(outdir,"/Filtering/Stats_tmp/Chunk_", i), full.names = T, pattern = "*_filtered_reads_automated.rds")
    filter_list[[i]]<-list()
    for (j in 1:length(filter_files)){
      filter_list[[i]][[j]]<-readRDS(filter_files[j])
    }
  }
  
  filter_list2<-list()
  for (i in 1:chunks){
    filter_list2[[i]]<-unlist(filter_list[i])
  }
  
  filter_comb<- unlist(filter_list2)
  
  write.table(filter_comb, paste0(outdir,"/Filtering/Retained_readIDs_automated.txt"), quote = F, sep="\n", row.names = F, col.names = F)
}

if (isTRUE(ncells_sel)){  
  filter_list<-list()
  for (i in 1:chunks){
    filter_files<-list.files(paste0(outdir,"/Filtering/Stats_tmp/Chunk_", i), full.names = T, pattern = "*cells.rds")
    filter_list[[i]]<-list()
    for (j in 1:length(filter_files)){
      filter_list[[i]][[j]]<-readRDS(filter_files[j])
    }
  }
  
  
  filter_list2<-list()
  for (i in 1:chunks){
    filter_list2[[i]]<-unlist(filter_list[i])
  }
  
  filter_comb<- unlist(filter_list2)
  
  write.table(filter_comb, paste0(outdir,"/Filtering/Retained_readIDs_",ncells ,"cells.txt"), quote = F, sep="\n", row.names = F, col.names = F)
}

if (isTRUE(custom_sel)){  
  filter_list<-list()
  for (i in 1:chunks){
    filter_files<-list.files(paste0(outdir,"/Filtering/Stats_tmp/Chunk_", i), full.names = T, pattern = "*_filtered_reads_custom.rds")
    filter_list[[i]]<-list()
    for (j in 1:length(filter_files)){
      filter_list[[i]][[j]]<-readRDS(filter_files[j])
    }
  }
  
  
  filter_list2<-list()
  for (i in 1:chunks){
    filter_list2[[i]]<-unlist(filter_list[i])
  }
  
  filter_comb<- unlist(filter_list2)
  
  write.table(filter_comb, paste0(outdir,"/Filtering/Retained_readIDs_custom.txt"), quote = F, sep="\n", row.names = F, col.names = F)
}
