#!/bin/bash
# IRD, Mivegec
# Main srunning script to run argentag pipeline version 2. Relying on taggy demux v 3.2.3
# Authors: Johannes W. Bagnoli
# Contact: johannes.bagnoli@ird.fr

# Params
yaml=$1






#sampleinfo:
sample=$(grep 'sample_name:' ${yaml} | awk '{print $2}')
input_path=$(grep 'input_path' ${yaml} | awk '{print $2}')

#general:
outdir=$(grep 'outdir' ${yaml} | awk '{print $2}')
nthreads=$(grep 'nthreads' ${yaml} | awk '{print $2}')
Start_stage=$(grep 'Start_stage' ${yaml} | awk '{print $2}')
End_stage=$(grep 'End_stage' ${yaml} | awk '{print $2}')

#taggy:
taggy_split=$(grep 'taggy_split:' ${yaml} | awk '{print $2}')
preset=$(grep 'preset' ${yaml} | awk '{print $2}')
taggy_params=$(grep 'taggy_params' ${yaml} | awk '{$1=""; print $0}')
keep_demux=$(grep 'keep_demux' ${yaml} | awk '{print $2}')

#filtering:
#Cellselection=$(grep 'Cellselection' ${yaml} | awk '{print $2}')
automated_sel=$(grep 'automated_sel' ${yaml} | awk '{print $2}')
ncells_sel=$(grep 'ncells_sel' ${yaml} | awk '{print $2}')
custom_sel=$(grep 'custom_sel' ${yaml} | awk '{print $2}')
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

## Checks
### Check if Rscript is excutable and Argentag pipline and yaml exists

if [ ! -f ${yaml} ]; then
  echo "Path to yaml file incorrect"
  exit 1
fi

if [[ ${yaml} != *.yaml ]] ; then
      echo "Path to yaml not a yaml file."
      exit 1
fi


if [ ! -d ${ArgenTAG_pipeline} ]; then
  echo "Path to ArgenTAG pipeline incorrect"
  exit 1
fi

if [ ! -x ${Rscript} ]; then
  echo "Rscript path is not executable"
  exit 1
fi

echo -e "\n....................Checking YAML file....................\n"
check=$(${Rscript} ${ArgenTAG_pipeline}/SCRIPTS/yaml_check.R ${yaml})
echo -e ${check}

if [[ ${check} == *"Final outcome..................PASS"* ]]; then
  echo "Pipeline is starting."
elif [[ ${check} == *"Final outcome..................FAIL"* ]]; then
  echo "Pipeline will not start."
  exit 1
else 
 echo "Yaml check not uscessful. Pipeline will not start."
 exit 1
fi

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


## Preprocessing
if [[ ${Start_ID} == 1 ]] ; then
  
  #create main output folder if it didn't exist
  if [[ ! -d ${outdir} ]] ; then
    mkdir -p ${outdir}
    if [ $? -ne 0 ] ; then
        echo "Please provide a valid output directory path."
        exit 1
    fi
  fi
  cd ${outdir}
  ## Unzip input fastq if ends in .gz
  if [[ ${input_path} == *.gz ]] ; then
    echo "unzipping input file"
    gzip -dkc ${input_path} > ${outdir}/${sample}.fastq
  elif [[ ${input_path} == *.fastq ]] ; then
    echo "copying input file"
    cp ${input_path}  ${outdir}/${sample}.fastq
  elif [[ ${input_path} == *.bam ]] ; then
    echo "Converting input file"
    ${samtoolsexc} fastq -@ ${nthreads} ${input_path} > ${outdir}/${sample}.fastq
  fi
  
  
  # Running taggy demux
  echo -e "\n....................Running taggy_demux....................\n"
  mkdir ${outdir}/taggy_demux
  if [[ ${taggy_split} > 1 ]] ; then
    mkdir ${outdir}/taggy_demux/split
    seqkit split2 ${outdir}/${sample}.fastq  -p ${taggy_split} -f -O ${outdir}/taggy_demux/split
    for (( i=1; i<=$taggy_split; i++ )); do mkdir ${outdir}/taggy_demux/split/taggy_demux_${i}  ; done
    cd ${ArgenTAG_pipeline}/taggy_demux_3.2.5d/taggy_demux-main/
    if [[ ${preset} == "none" && ${taggy_params} == "none" ]] ; then
      for (( i=1; i<=$taggy_split; i++ )); do bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/split/taggy_demux_${i} -s ${outdir}/taggy_demux/split/${sample}.part_00${i}.fastq ; done
    elif [[ ${preset} != "none" && ${taggy_params} == "none" ]] ; then
      for (( i=1; i<=$taggy_split; i++ )); do bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/split/taggy_demux_${i}  -s ${outdir}/taggy_demux/split/${sample}.part_00${i}.fastq --presets=${preset} ; done
    elif [[ ${preset} == "none" && ${taggy_params} != "none" ]] ; then
      for (( i=1; i<=$taggy_split; i++ )); do bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/split/taggy_demux_${i}  -s ${outdir}/taggy_demux/split/${sample}.part_00${i}.fastq ${taggy_params} ; done
    else
      for (( i=1; i<=$taggy_split; i++ )); do bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/split/taggy_demux_${i}  -s ${outdir}/taggy_demux/split/${sample}.part_00${i}.fastq --presets=${preset} ${taggy_params} ; done
    fi
    rm ${outdir}/taggy_demux/split/*.fastq
    mkdir  ${outdir}/taggy_demux/fastq
    cd ${outdir}/taggy_demux/split/taggy_demux_1/fastq/
    for i in  *".fastq" ; do cat "${outdir}/taggy_demux/split/taggy_demux_"*"/fastq/$i" > "${outdir}/taggy_demux/fastq/$i" ; done
    cd ${outdir}
    for (( i=1; i<=$taggy_split; i++ )); do rm -r ${outdir}/taggy_demux/split/taggy_demux_${i}/fastq  ; done
  else
    cd ${ArgenTAG_pipeline}/taggy_demux_3.2.5d/taggy_demux-main/
    if [[ ${preset} == "none" && ${taggy_params} == "none" ]] ; then
      bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq 
    elif [[ ${preset} != "none" && ${taggy_params} == "none" ]] ; then
      bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq --presets=${preset}
    elif [[ ${preset} == "none" && ${taggy_params} != "none" ]] ; then
      bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq ${taggy_params}
    else
      bin/taggy_demux -T ${nthreads} -o ${outdir}/taggy_demux/ -s ${outdir}/${sample}.fastq --presets=${preset} ${taggy_params}
    fi
  fi
  rm ${outdir}/${sample}.fastq
  pigz ${outdir}/taggy_demux/fastq/*fastq
  cat ${outdir}/taggy_demux/fastq/*fastq.gz  > ${outdir}/taggy_demux/fastq/${sample}.taggydemux.fastq.gz
fi

# Stats
if [[ ${Start_ID} < 3 && ${End_ID} > 1 ]] ; then
  echo -e "\n....................Running Cell Selection....................\n"
  ${Rscript} ${ArgenTAG_pipeline}/SCRIPTS/Stats_v2.R ${yaml}
  
  if [[ ${keep_temp_stats} == "no" ]] ; then
    rm -r ${outdir}/Filtering/Stats_tmp
  fi
  
  if [[ ${keep_demux} == "no" ]] ; then
    find ${outdir}/taggy_demux/fastq/ -type f -name "[0-9][0-9][0-9][0-9].fastq.gz" -delete
  fi
fi



# Mapping
if [[ ${Start_ID} < 4 && ${End_ID} > 2 ]] ; then
  echo -e "\n....................Running Mapping....................\n"
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
  echo -e "\n....................Running Filtering....................\n"
  if [[ ${automated_sel} ==  "yes" ]] ; then
    ${samtoolsexc} view -@ ${nthreads} -b --qname-file ${outdir}/Filtering/Retained_readIDs_automated.txt ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.bam  > ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam 
    ${samtoolsexc} index -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam
    ${samtoolsexc} flagstats -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam > ${outdir}/Mapping/flagstats/${sample}.taggydemux.mapped.${genome_name}.filt_auto.flagstats.txt 
  fi
  
  if [[ ${ncells_sel} ==  "yes" ]] ; then
    echo "${ncells} targeted cells in yaml file"
    ncells=$(cat  ${outdir}/Filtering/updated_ncells.txt)
    echo "${ncells} cells will be used for Filtering"
    ${samtoolsexc} view -@ ${nthreads} -b --qname-file ${outdir}/Filtering/Retained_readIDs_"${ncells}"cells.txt ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.bam  > ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam 
    ${samtoolsexc} index -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam 
    ${samtoolsexc} flagstats -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam  > ${outdir}/Mapping/flagstats/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.flagstats.txt 
  fi
  
  if [[ ${custom_sel} ==  "yes" ]] ; then
    ${samtoolsexc} view -@ ${nthreads} -b --qname-file ${outdir}/Filtering/Retained_readIDs_custom.txt ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.bam  > ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_custom.bam 
    ${samtoolsexc} index -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_custom.bam
    ${samtoolsexc} flagstats -@ ${nthreads} ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_custom.bam > ${outdir}/Mapping/flagstats/${sample}.taggydemux.mapped.${genome_name}.filt_custom.flagstats.txt 
  fi
fi

# Bambu
if [[ ${Start_ID} < 6 && ${End_ID} > 4 ]] ; then
  if [[ ${automated_sel} ==  "yes" ]] ; then
    echo -e "\n....................Running Bambu for automated cell selection....................\n"
    ${Rscript}  ${ArgenTAG_pipeline}/SCRIPTS/bambu_v2.R ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_auto.bam  ${yaml} automated
  fi
  
  if [[ ${ncells_sel} ==  "yes" ]] ; then
    ncells=$(cat  ${outdir}/Filtering/updated_ncells.txt)
    echo -e "\n....................Running Bambu for ${ncells} cells selection....................\n"
    ${Rscript}  ${ArgenTAG_pipeline}/SCRIPTS/bambu_v2.R ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_"${ncells}"cells.bam  ${yaml} "${ncells}cells"
  fi
  
  if [[ ${custom_sel} ==  "yes" ]] ; then
    echo -e "\n....................Running Bambu for custom cell selection....................\n"
    ${Rscript}  ${ArgenTAG_pipeline}/SCRIPTS/bambu_v2.R ${outdir}/Mapping/${sample}.taggydemux.mapped.${genome_name}.sorted.filt_custom.bam  ${yaml} custom
  fi
fi


# END
echo "DONE"

