**yaml paramters:** <br />
--- <br />
**sampleinfo:** <br />
  sample: character, name of the sample <br />
  input_path: path, input file, fastq, fastq.gz or bam <br />
**general:** <br />
  outdir: path, output directory, if does not exist will create <br />
  nthreads: nr of threads to use <br />
  nchunks: nr of chunks to use during stats, more chunks lead to higher RAM usage, 15 max recommended <br />
  Start_stage: which part to start with, one of taggy, stats,  mapping, filtering, bambu <br />
  End_stage: last parts to run included ,bambu taggy, stats,  mapping, filtering, bambu <br />
**taggy:** <br />
  taggy_split: number of seperate taggy_demux runs, numeric 1-9  <br />
  preset: taggy demux preset paramater, e.g. ont+v1 <br />
  keep_demux: TRUE/FALSE keep seperate taggy_demux fastq files <br />
**filtering:** <br />
  Cellselection: auto, both or ncells, if auto, ncells parameter will be ignored. <br />
  ncells: number of target cells, e.g.10000, do not set higher than 30000, can lead to excessive RAM usage during Bambu <br />
  minreadlength: minimal mean read length for a barcode to be considered, normally 400-600 <br />
  keep_temp_stats: keep temporary stats files <br />
**reference:** <br />
  genome_name: name of the genome, will be added to mapped files to distiungish them <br />
  genome_path: path to genome fasta <br />
  gtf_path: path to genome gtf file  <br />
  bed_path: path to genome bed file for minimap <br />
  splice_size: minimap maximum splice size <br />
  gene_keys: Identifying character for genes of species, e.g. ENSG for human, if several genomes were used seperate by , e.g. ENGS, ENSMUS <br />
  species_names: names of species, e.g.  human, if several genomes were used seperate by , e.g. human , mouse <br /> <br />
**tools:** paths for tools used <br />
  Rscript: Rscript <br />
  ArgenTAG_pipeline: /home/Argentag_Pipeline_V2 <br />
  samtoolsexc: samtools <br />
  minimap2_exc: minimap2 <br />
--- <br />
