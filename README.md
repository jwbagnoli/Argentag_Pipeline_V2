yaml paramters:
---
sampleinfo: 
  sample: character, name of the sample
  input_path: path, input file, fastq, fastq.gz or bam
general:
  outdir: path, output directory, if does not exist will create
  nthreads: nr of threads to use
  nchunks: nr of chunks to use during stats, more chunks lead to higher RAM usage, 15 max recommended
  Start_stage: which part to start with, one of taggy, stats,  mapping, filtering, bambu
  End_stage: last parts to run included ,bambu taggy, stats,  mapping, filtering, bambu
taggy:
  taggy_split: number of seperate taggy_demux runs, numric 1-9, 
  preset: taggy demux preset paramater, e.g. ont+v1
  keep_demux: TRUE/FALSE keep seperate taggy_demux fastq files
filtering:
  Cellselection: auto, both or ncells, if auto, ncells parameter will be ignored.
  ncells: number of target cells, e.g.10000, do not set higher than 30000, can lead to excessive RAM usage during Bambu
  minreadlength: minimal mean read length for a barcode to be considered, normally 400-600
reference:
  genome_name: name of the genome, will be added to mapped files to distiungish them
  genome_path: path to genome fasta
  gtf_path: path to genome gtf file 
  bed_path: path to genome bed file for minimap
  splice_size: minimap maximum splice size
tools: paths to tools used
  Rscript: /usr/bin/Rscript
  ArgenTAG_pipeline: /home/Tools/Argentag_Pipeline_V2
  samtoolsexc: /home/Tools/samtools-1.23/samtools
  minimap2_exc: minimap2
---
