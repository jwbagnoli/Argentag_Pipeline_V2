**Single sample:** <br />
<<<<<<< HEAD
Usage: bash Argentag_Pipeline_V2.sh path/singlesample.yaml
=======
Usage: bash Argentag_Pipeline_V2.sh path/single.yaml
>>>>>>> 7faa725278ce442fa681cac37ed5e3c2ba246be4

Make sure to use singlesample.yaml template!!!

**yaml paramters:** <br />
--- <br />
**multi:** no <br />
**sampleinfo:** <br />
  sample_name: character, name of the sample <br />
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
  keep_demux: yes/no keep seperate taggy_demux fastq files <br />
**filtering:** several cell selections can be run. mapping is performed on unfiltered reads.<br />
  automated_sel: yes/no, run automated cell selection, ignores ncells, minreads and minreadlength <br />
  ncells_sel: yes/no, run ncells cell selection, ignores minreads. filters for minreadlength and selects  ncells barcodes with highest number of reads <br />
  custom_sel: yes/no, runs custom cell selection, ignores ncells, selects barcodes >= minreads and >=  minreadlength <br />
  ncells: number of target cells, e.g.10000, do not set higher than 30000, can lead to excessive RAM usage during Bambu, ignored in automated and custom cell selection <br />
  minreads: minimal number of reads for a barcode to be considered, normally 400-1000, ignored in automated and ncell cell selection <br />
  minreadlength: minimal mean read length for a barcode to be considered, normally 400-600, ignored in automated cell selection <br />
  keep_temp_stats: yes/no keep temporary stats files <br />
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


*Multiple samples:** <br />
<<<<<<< HEAD
Usage: bash Argentag_Pipeline_V2.sh path/multisample.yaml

Make sure to use multisample.yaml template!!! 
Each sample needs to have it's own subsection for filtering and reference sections. 
Each parameter is completely independent for each sample making it possible to run samples from different species, or wastly different sequencing qualities.
The sampleinfo section is structure differently as well!
The sample specific folders and newly created yaml files will be created automatically.
Example shows 2 samples!
=======
Usage: bash Argentag_Pipeline_V2.sh path/multi.yaml

Make sure to use multisample.yaml template!!! Each sample needs to have it's own subsection for filtering and reference sections. the sampleinfo section is structure differently as well!
Example show 2 samples! 
>>>>>>> 7faa725278ce442fa681cac37ed5e3c2ba246be4

**yaml paramters:** <br />
--- <br />
**multi:** yes <br />
**sampleinfo:** <br />
  mysample1: path, input file for sample 1, fastq, fastq.gz or bam <br />
  mysample2: path, input file for sample 2, fastq, fastq.gz or bam <br />
**general:** <br />
  outdir: path, output directory, if does not exist will create <br />
  nthreads: nr of threads to use <br />
  nchunks: nr of chunks to use during stats, more chunks lead to higher RAM usage, 15 max recommended <br />
  Start_stage: which part to start with, one of taggy, stats,  mapping, filtering, bambu <br />
  End_stage: last parts to run included ,bambu taggy, stats,  mapping, filtering, bambu <br />
**taggy:** <br />
  taggy_split: number of seperate taggy_demux runs, numeric 1-9  <br />
  preset: taggy demux preset paramater, e.g. ont+v1 <br />
  keep_demux: yes/no keep seperate taggy_demux fastq files <br />
**filtering:** several cell selections can be run. mapping is performed on unfiltered reads.<br />
  mysample1: <br />
    automated_sel: yes/no, run automated cell selection, ignores ncells, minreads and minreadlength <br />
    ncells_sel: yes/no, run ncells cell selection, ignores minreads. filters for minreadlength and selects  ncells barcodes with highest number of reads <br />
    custom_sel: yes/no, runs custom cell selection, ignores ncells, selects barcodes >= minreads and >=  minreadlength <br />
    ncells: number of target cells, e.g.10000, do not set higher than 30000, can lead to excessive RAM usage during Bambu, ignored in automated and custom cell selection <br />
    minreads: minimal number of reads for a barcode to be considered, normally 400-1000, ignored in automated and ncell cell selection <br />
    minreadlength: minimal mean read length for a barcode to be considered, normally 400-600, ignored in automated cell selection <br />
    keep_temp_stats: yes/no keep temporary stats files <br />
  mysample2: <br />
    automated_sel: yes/no, run automated cell selection, ignores ncells, minreads and minreadlength <br />
    ncells_sel: yes/no, run ncells cell selection, ignores minreads. filters for minreadlength and selects  ncells barcodes with highest number of reads <br />
    custom_sel: yes/no, runs custom cell selection, ignores ncells, selects barcodes >= minreads and >=  minreadlength <br />
    ncells: number of target cells, e.g.10000, do not set higher than 30000, can lead to excessive RAM usage during Bambu, ignored in automated and custom cell selection <br />
    minreads: minimal number of reads for a barcode to be considered, normally 400-1000, ignored in automated and ncell cell selection <br />
    minreadlength: minimal mean read length for a barcode to be considered, normally 400-600, ignored in automated cell selection <br />
    keep_temp_stats: yes/no keep temporary stats files <br />
**reference:** <br />
  mysample1: <br />
    genome_name: name of the genome, will be added to mapped files to distiungish them <br />
    genome_path: path to genome fasta <br />
    gtf_path: path to genome gtf file  <br />
    bed_path: path to genome bed file for minimap <br />
    splice_size: minimap maximum splice size <br />
    gene_keys: Identifying character for genes of species, e.g. ENSG for human, if several genomes were used seperate by , e.g. ENGS, ENSMUS <br />
    species_names: names of species, e.g.  human, if several genomes were used seperate by , e.g. human , mouse <br /> <br />
  mysample2: <br />
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