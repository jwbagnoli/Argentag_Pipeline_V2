#!/bin/bash
# IRD, Mivegec
# Main srunning script to run argentag pipeline version 2. Relying on taggy demux v 3.2.3
# Authors: Johannes W. Bagnoli
# Contact: johannes.bagnoli@ird.fr

# Params
yaml=$1

sample=$(grep 'sample:' ${yaml} | awk '{print $2}')
fastqc_path=$(grep 'fastqc_path' ${yaml} | awk '{print $2}')
nthreads=$(grep 'nthreads' ${yaml} | awk '{print $2}')
outdir=$(grep 'outdir' ${yaml} | awk '{print $2}')
preset=$(grep 'split' ${yaml} | awk '{print $2}')
genome_path=$(grep 'genome_path' ${yaml} | awk '{print $2}')
gtf_path=$(grep 'gtf_path' ${yaml} | awk '{print $2}')
bed_path=$(grep 'bed_path' ${yaml} | awk '{print $2}')
genome_name=$(grep 'genome_name' ${yaml} | awk '{print $2}')
splice_size=$(grep 'splice_size:' ${yaml} | awk '{print $2}')
Cellselection=$(grep 'Cellselection' ${yaml} | awk '{print $2}')
ncells=$(grep 'ncells:' ${yaml} | awk '{print $2}')
Rscript=$(grep 'Rscript' ${yaml} | awk '{print $2}')
ArgenTAG_pipeline=$(grep 'ArgenTAG_pipeline' ${yaml} | awk '{print $2}')
samtoolsexc=$(grep 'samtoolsexc' ${yaml} | awk '{print $2}')
minimap2_exc=$(grep 'minimap2_exc' ${yaml} | awk '{print $2}')

## Preprocessing
#create main output folder if it didn't exist
if [[ ! -d ${outdir} ]] ; then
  mkdir -p ${outdir}
  if [ $? -ne 0 ] ; then
      echo "Please provide a valid output directory path."
      exit 1
  fi
fi

echo "unzipping input file"
gzip -dkc ${fastqc_path} > ${outdir}/${sample}.fastq

# Running taggy demux
echo "Running taggy_demux"
mkdir ${outdir}/taggy_demux
cd ${ArgenTAG_pipeline}/taggy_demux_3.2.3/taggy_demux-main/
bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq --presets=${preset}

rm ${outdir}/${sample}.fastq
pigz ${outdir}/taggy_demux/fastq/*fastq
cat ${outdir}/taggy_demux/fastq/*fastq.gz  > ${outdir}/taggy_demux/fastq/${sample}.taggydemux.fastq.gz

# Stats
echo "Running Cell Selection"
${Rscript} ${ArgenTAG_pipeline}/SCRIPTS/Stats_v2.R ${yaml}
rm -r ${outdir}/Filtering/Stats_tmp
find ${outdir}/taggy_demux/fastq/ -type f -name "[0-9][0-9][0-9][0-9].fastq.gz" -delete

# Mapping
echo "Running Mapping"
mkdir ${outdir}/Mapping
${minimap2_exc} -t ${nthreads} -ax splice -G ${splice_size} --junc-bed ${bed_path} ${genome_path} ${outdir}/taggy_demux/fastq/${sample}.taggydemux.fastq.gz > ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sam
${samtoolsexc} view -@ ${nthreads} -b -S ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sam > ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.bam
rm ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sam


mkdir ${outdir}/Mapping/flagstats
${samtoolsexc} flagstats -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.bam > ${outdir}/Mapping/flagstats/${sample}.taggydemux.mapped.${genome_name}.flagstats.txt 

${samtoolsexc} sort -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.bam > ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.bam
${samtoolsexc} index -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.bam
rm ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.bam 

# Filtering Mapping
echo "Running Filtering"
if [[ ${Cellselection} ==  "auto" || ${Cellselection} ==  "both" ]] ; then
  ${samtoolsexc} view -@ ${nthreads} -b --qname-file ${outdir}/Filtering/Retained_readIDs_automated.txt ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.bam  > ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam 
  ${samtoolsexc} index -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam
  ${samtoolsexc} flagstats -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam > ${outdir}/Mapping/flagstats/${sample}.taggydemux.mapped.${genome_name}.filt_auto.flagstats.txt 
fi

if [[ ${Cellselection} ==  "ncells" || ${Cellselection} ==  "both" ]] ; then
  ${samtoolsexc} view -@ ${nthreads} -b --qname-file ${outdir}/Filtering/Retained_readIDs_"${ncells}"cells.txt ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.bam  > ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam 
  ${samtoolsexc} index -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam 
  ${samtoolsexc} flagstats -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam  > ${outdir}/Mapping/flagstats/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.flagstats.txt 
fi


# Bambu
if [[ ${Cellselection} ==  "auto" || ${Cellselection} ==  "both" ]] ; then
  echo "Running Bambu for automated cell selection"
  ${Rscript}  ${ArgenTAG_pipeline}/SCRIPTS/bambu_v2.R ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam  ${yaml} automated
fi

if [[ ${Cellselection} ==  "ncells" || ${Cellselection} ==  "both" ]] ; then
  echo "Running Bambu for ${ncells} cells selection"
  ${Rscript}  ${ArgenTAG_pipeline}/SCRIPTS/bambu_v2.R ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam  ${yaml} "${ncells}cells"
fi


# END
echo "DONE"

