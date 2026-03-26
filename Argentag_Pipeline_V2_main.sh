#!/bin/bash
# IRD, Mivegec
# Main srunning script to run argentag pipeline version 2. Relying on taggy demux v 3.2.3
# Authors: Johannes W. Bagnoli
# Contact: johannes.bagnoli@ird.fr

# Params
yaml=$1

#sampleinfo:
sample=$(grep 'sample:' ${yaml} | awk '{print $2}')
fastq_path=$(grep 'fastq_path' ${yaml} | awk '{print $2}')

#general:
outdir=$(grep 'outdir' ${yaml} | awk '{print $2}')
nthreads=$(grep 'nthreads' ${yaml} | awk '{print $2}')
Start_stage=$(grep 'Start_stage' ${yaml} | awk '{print $2}')
End_stage=$(grep 'End_stage' ${yaml} | awk '{print $2}')

#taggy:
preset=$(grep 'preset' ${yaml} | awk '{print $2}')
taggy_params=$(grep 'taggy_params' ${yaml} | awk '{print $2}')
keep_demux=$(grep 'keep_demux' ${yaml} | awk '{print $2}')

#filtering:
Cellselection=$(grep 'Cellselection' ${yaml} | awk '{print $2}')
ncells=$(grep 'ncells:' ${yaml} | awk '{print $2}')
keep_temp_stats=$(grep 'keep_temp_stats' ${yaml} | awk '{print $2}')

#reference:
genome_name=$(grep 'genome_name' ${yaml} | awk '{print $2}')
genome_path=$(grep 'genome_path' ${yaml} | awk '{print $2}')
gtf_path=$(grep 'gtf_path' ${yaml} | awk '{print $2}')
bed_path=$(grep 'bed_path' ${yaml} | awk '{print $2}')
splice_size=$(grep 'splice_size:' ${yaml} | awk '{print $2}')

#tools
Rscript=$(grep 'Rscript' ${yaml} | awk '{print $2}')
ArgenTAG_pipeline=$(grep 'ArgenTAG_pipeline' ${yaml} | awk '{print $2}')
samtoolsexc=$(grep 'samtoolsexc' ${yaml} | awk '{print $2}')
minimap2_exc=$(grep 'minimap2_exc' ${yaml} | awk '{print $2}')



## Encoding stages

if [[ ${Start_stage}  == "taggy" ]] ; then
      Start_ID=1
elif [[ ${Start_stage}  == "stats" ]] ; then
    Start_ID=2
elif [[ ${Start_stage}  == "mapping" ]] ; then
    Start_ID=3
elif [[ ${Start_stage}  == "filtering" ]] ; then
    Start_ID=4
elif [[ ${Start_stage}  == "bambu" ]] ; then
    Start_ID=5
fi


if [[ ${End_stage}  == "taggy" ]] ; then
      End_ID=1
elif [[ ${End_stage}  == "stats" ]] ; then
    End_ID=2
elif [[ ${End_stage}  == "mapping" ]] ; then
    End_ID=3
elif [[ ${End_stage}  == "filtering" ]] ; then
    End_ID=4
elif [[ ${End_stage}  == "bambu" ]] ; then
    End_ID=5
fi

## Checks
### Check if input file is fastq or fastq.gz, if not exit
if [[ ${fastq_path} != *fastq.gz && ${fastq_path} != *.fastq ]] ; then
      echo "Please provide a fastq or fastq.gz input file."
      exit 1
fi

if [[ ${End_ID} <  ${Start_ID} ]] ; then
      echo "End stage comes before Start stage, please provide and End Stage which is downstream of the Start stage."
      exit 1
fi


## Preprocessing
#create main output folder if it didn't exist

if [[ ${Start_ID} == 1 ]] ; then
  if [[ ! -d ${outdir} ]] ; then
    mkdir -p ${outdir}
    if [ $? -ne 0 ] ; then
        echo "Please provide a valid output directory path."
        exit 1
    fi
  fi

  ## Unzip input fastq if ends in .gz
  if [[ ${fastq_path} == *.gz ]] ; then
    echo "unzipping input file"
    gzip -dkc ${fastq_path} > ${outdir}/${sample}.fastq
  else
    echo "copying input file"
    cp ${fastq_path}  ${outdir}/${sample}.fastq
  fi
  
  
  # Running taggy demux
  echo "Running taggy_demux"
  mkdir ${outdir}/taggy_demux
  cd ${ArgenTAG_pipeline}/taggy_demux_3.2.3/taggy_demux-main/
  if [[ ${preset} == "none" && ${taggy_params} == "none" ]] ; then
    bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq 
  elif [[ ${preset} != "none" && ${taggy_params} == "none" ]] ; then
    bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq --presets=${preset}
  elif [[ ${preset} == "none" && ${taggy_params} != "none" ]] ; then
    bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq ${taggy_params}
  else
    bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq --presets=${preset} ${taggy_params}
  fi

  rm ${outdir}/${sample}.fastq
  pigz ${outdir}/taggy_demux/fastq/*fastq
  cat ${outdir}/taggy_demux/fastq/*fastq.gz  > ${outdir}/taggy_demux/fastq/${sample}.taggydemux.fastq.gz
fi


# Stats
if [[ ${Start_ID} < 3 && ${End_ID} > 1 ]] ; then
  echo "Running Cell Selection"
  ${Rscript} ${ArgenTAG_pipeline}/SCRIPTS/Stats_v2.R ${yaml}
  
  if [[ ${keep_temp_stats} == FALSE ]] ; then
    rm -r ${outdir}/Filtering/Stats_tmpe
  fi
  
  if [[ ${keep_demux} == FALSE ]] ; then
    find ${outdir}/taggy_demux/fastq/ -type f -name "[0-9][0-9][0-9][0-9].fastq.gz" -delete
  fi
fi


# Mapping
if [[ ${Start_ID} < 4 && ${End_ID} > 2 ]] ; then
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
fi

# Filtering Mapping
if [[ ${Start_ID} < 5 && ${End_ID} > 3 ]] ; then
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
fi

# Bambu
if [[ ${Start_ID} < 6 && ${End_ID} > 4 ]] ; then
  if [[ ${Cellselection} ==  "auto" || ${Cellselection} ==  "both" ]] ; then
    echo "Running Bambu for automated cell selection"
    ${Rscript}  ${ArgenTAG_pipeline}/SCRIPTS/bambu_v2.R ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam  ${yaml} automated
  fi
  
  if [[ ${Cellselection} ==  "ncells" || ${Cellselection} ==  "both" ]] ; then
    echo "Running Bambu for ${ncells} cells selection"
    ${Rscript}  ${ArgenTAG_pipeline}/SCRIPTS/bambu_v2.R ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam  ${yaml} "${ncells}cells"
  fi
fi

# END
echo "DONE"

