#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


suppressPackageStartupMessages(library(dplyr, quietly = T, verbose = F))

yaml_path<- args[1]
#yaml_path<-"/home/TALMAN/Tools/Argentag_Pipeline_V2/test/test_errors.yaml"

ym<-yaml::read_yaml(yaml_path)




FAILS<-c()
errormsg<-c()
#### sample Infos ----

# Check if sample name is integer
if (!(is.character(ym$sampleinfo$sample))){
  errormsg_tmp<-"sample name is not a character"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[1]<-T
}else{
  FAILS[1]<-F
}



# Check if input file exists
if (sum(grep(pattern="*.fastq.gz$", ym$sampleinfo$input_path),
        grep(pattern="*.fastq$", ym$sampleinfo$input_path),
        grep(pattern="*.bam$", ym$sampleinfo$input_path)) == 0){
  errormsg_tmp<-"Input file not fastq, fastq.gz or bam"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[2]<-T
}else{
  FAILS[2]<-F 
}


if (!(file.exists(ym$sampleinfo$input_path))){
  errormsg_tmp<-"input file NOT found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[3]<-T
}else{
  FAILS[3]<-F 
  }


#### general ----
stages<- c(1,2,3,4,5)
names(stages)<- c("taggy", "stats",  "mapping", "filtering", "bambu")
# Check if output directory is possible (parent is checked)
out_parent<-paste(strsplit(ym$general$outdir, split = "/")[[1]][1:(length(strsplit(ym$general$outdir, split = "/")[[1]])-1)], collapse = "/")
if (!(file.exists(out_parent))){
  errormsg_tmp<-"output path NOT found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[4]<-T
}else{
  FAILS[4]<-F 
}

# Check if number of threads is integer and does not exceed system cores
if (!(is.numeric(ym$general$nthreads))){
  errormsg_tmp<-"Number of threads is not integer"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[5]<-T
}else if(is.numeric(ym$general$nthreads) & ym$general$nthreads > parallel::detectCores()){
  errormsg_tmp<-paste0("Number of threads is larger than avaiable cores (",parallel::detectCores(), ")")
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[5]<-T
}else{
  FAILS[5]<-F
}

# Check if number of chunks is integer
if (!(is.numeric(ym$general$nchunks))){
  errormsg_tmp<-"Number of chunks is not integer"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[6]<-T
}else{
  FAILS[6]<-F
}

# Check if Start stage is possible
if (!(ym$general$Start_stage %in% names(stages))){
  errormsg_tmp<-paste0("Start stage is not one of ", paste(names(stages), collapse = ", "))
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[7]<-T
}else{
  FAILS[7]<-F
}

# Check if End stage is possible and not before Start stage 
if (!(ym$general$End_stage %in% names(stages))){
  errormsg_tmp<-paste0("End stage is not one of ",  paste(names(stages), collapse = ", "))
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[8]<-T
}else if(ym$general$End_stage %in% names(stages) & as.numeric(stages[ym$general$End_stage]) < as.numeric(stages[ym$general$Start_stage]) ){
  errormsg_tmp<-"End stage is before than Start Stage"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[8]<-T
}else{
  FAILS[8]<-F
}


#### taggy ----
# Check if splitting of taggy demux is integer and maximum 9
if (!(is.numeric(ym$taggy$taggy_split))){
  errormsg_tmp<-"Number of splitting for tagging demux is not integer"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[9]<-T
}else if(is.numeric(ym$taggy$taggy_split) & ym$taggy$taggy_split > 9){
  errormsg_tmp<-"Number of splitting bigger than 9"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[9]<-T
}else{
  FAILS[9]<-F
}

# Check if preset is ok
versions<-c("v1", "v2")
plats<-c("hifi", "ont", "illu")
presets<-c(as.character(sapply(versions, FUN = function(x) paste0(x,"+", plats))), as.character(sapply(plats, FUN = function(x) paste0(x,"+", versions))))
preset_print<-paste(presets, collapse = ",")

if (!(ym$taggy$preset %in% presets)){
  errormsg_tmp<-paste0(".....preset for tagging demux is not compatible with options (",preset_print,")")
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[10]<-T
}else{
  FAILS[10]<-F
}

# Check if keep demux is logical
if (!(is.logical(ym$taggy$keep_demux))){
  errormsg_tmp<-"Keeping demultiplexed files must be logical (TRUE/FALSE)"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[11]<-T
}else{
  FAILS[11]<-F
}


#### filtering ----
# Check if automated cell selection is logical
if (!(is.logical(ym$filtering$automated_sel))){
  errormsg_tmp<-"Automated cell selection (automated_sel) must be logical (TRUE/FALSE)"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[12]<-T
}else{
  FAILS[12]<-F
}

# Check if ncells cell selection is logical
if (!(is.logical(ym$filtering$ncells_sel))){
  errormsg_tmp<-"Number of cells cell selection (ncells_sel) must be logical (TRUE/FALSE)"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[13]<-T
}else{
  FAILS[13]<-F
}

# Check if custom cell selection is logical
if (!(is.logical(ym$filtering$custom_sel))){
  errormsg_tmp<-"Custom cell selection (custom_sel) must be logical (TRUE/FALSE)"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[14]<-T
}else{
  FAILS[14]<-F
}

# Check if ncells is numerical
if (!(is.numeric(ym$filtering$ncells))){
  errormsg_tmp<-"Number of cells (ncells) must be numerical"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[15]<-T
}else{
  FAILS[15]<-F
}

# Check if minreads is numerical
if (!(is.numeric(ym$filtering$minreads))){
  errormsg_tmp<-"Minimal number of reads (minreads) must be numerical"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[16]<-T
}else{
  FAILS[16]<-F
}

# Check if minreadlength is numerical
if (!(is.numeric(ym$filtering$minreadlength))){
  errormsg_tmp<-"Minimal read length (minreadlength) must be numerical"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[17]<-T
}else{
  FAILS[17]<-F
}

# Check if keep temporary stats is logical
if (!(is.logical(ym$filtering$keep_temp_stats))){
  errormsg_tmp<-"Keeping temporary stats files must be logical must be logical (TRUE/FALSE)"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[18]<-T
}else{
  FAILS[18]<-F
}

#### Reference ----
# Check if genome name is character
if (!(is.character(ym$reference$genome_name))){
  errormsg_tmp<-"Genome name is not a character"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[19]<-T
}else{
  FAILS[19]<-F
}

# Check if genome path exists
if (!(file.exists(ym$reference$genome_path))){
  errormsg_tmp<-"Genome path not found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[20]<-T
}else{
  FAILS[20]<-F
}

# Check if gtf path exists
if (!(file.exists(ym$reference$gtf_path))){
  errormsg_tmp<-"GTF annotation path not found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[21]<-T
}else{
  FAILS[21]<-F
}

# Check if bed path exists
if (!(file.exists(ym$reference$bed_path))){
  errormsg_tmp<-"BED annotation path not found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[22]<-T
}else{
  FAILS[22]<-F
}

# Check if splice size is numeric
if (!(is.numeric(ym$reference$splice_size))){
  errormsg_tmp<-"Splice size is not numeric"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[23]<-T
}else{
  FAILS[23]<-F
}

# Check if gene keys is character
if (!(is.character(ym$reference$gene_keys))){
  errormsg_tmp<-"gene keys are not charcter"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[24]<-T
}else{
  FAILS[24]<-F
}

# Check if species name is character
if (!(is.character(ym$reference$species_names))){
  errormsg_tmp<-"species names are not a character"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[25]<-T
}else{
  FAILS[25]<-F
}

# Check if species name and gene keys are the same size
if (!(length(ym$reference$species_names) == length(ym$reference$gene_keys))){
  errormsg_tmp<-"species names and gene keys are not the same length"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[26]<-T
}else{
  FAILS[26]<-F
}

#### Tools ----
is_bin_on_path = function(bin) {
  exit_code = suppressWarnings(system2("command", args = c("-v", bin), stdout = FALSE))
  return(exit_code == 0)
}

# Check if Rscript exists
if (!(file.exists(ym$tools$Rscript)) & !(is_bin_on_path(ym$tools$Rscript))){
  errormsg_tmp<-"Rscript not found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[27]<-T
}else{
  FAILS[27]<-F
}

# Check if Argentag piepline exists
if (!(file.exists(ym$tools$ArgenTAG_pipeline))){
  errormsg_tmp<-"ArgenTAG pipeline not found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[28]<-T
}else{
  FAILS[28]<-F
}

# Check if samtools exists
if (!(file.exists(ym$tools$samtoolsexc)) & !(is_bin_on_path(ym$tools$samtoolsexc))){
  errormsg_tmp<-"samtools not found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[29]<-T
}else{
  FAILS[29]<-F
}

# Check if minimap2 exists
if (!(file.exists(ym$tools$minimap2_exc)) & !(is_bin_on_path(ym$tools$minimap2_exc))){
  errormsg_tmp<-"minimap2 not found"
  errormsg<-c(errormsg, errormsg_tmp)
  FAILS[30]<-T
}else{
  FAILS[30]<-F
}

#### Final count ----
failcounts<-sum(FAILS)

if (failcounts == 0){
  decision<-"PASS"
  out<-data.frame(description= c("Final outcome", "Nr.Errors"),
                  Value=c(decision, failcounts))
  paste0("\n",out$description[1], "..................", out$Value[1], "\n",
             out$description[2], "......................", out$Value[2], "\n\n",
             "No errors in yaml.\n")

}else{
  decision<-"FAIL"
  out<-data.frame(description= c("Final outcome", "Nr.Errors", "Errors", rep("", length(errormsg))),
                  Value=c(decision, failcounts, "",errormsg ))

  paste0("\n",out$description[1], "..................", out$Value[1], "\n",
             out$description[2], "......................", out$Value[2],"\n\n",
             out$description[3], ":\n",
             paste(errormsg, collapse = "\n"),
             "\n\n",failcounts, " error(s) found in yaml.\n")

}



