#!/bin/bash
# IRD, Mivegec
# Setup script to run argentag pipeline version 2.
# Authors: Johannes W. Bagnoli
# Contact: johannes.bagnoli@ird.fr

# Params
yaml=$1

# Multisample
multi=$(grep 'multi' ${yaml} | awk '{print $2}')

# tools
Rscript=$(grep 'Rscript' ${yaml} | awk '{print $2}')
ArgenTAG_pipeline=$(grep 'ArgenTAG_pipeline' ${yaml} | awk '{print $2}')

#general:
outdir=$(grep 'outdir' ${yaml} | awk '{print $2}')
Start_stage=$(grep 'Start_stage' ${yaml} | awk '{print $2}')

## Checks
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


# Run pipeline
if [[ ${multi} == "yes" ]] ; then
  echo -e "\n....................Preparing Multiple sample pipeline....................\n"
  if [[ ${Start_stage} == "taggy" ]] ; then
    if [[ ! -d ${outdir} ]] ; then
      mkdir -p ${outdir}
      if [ $? -ne 0 ] ; then
          echo "Please provide a valid output directory path."
          exit 1
      fi
    fi
  fi
  cd ${outdir}
  ${Rscript} ${ArgenTAG_pipeline}/SCRIPTS/MultiPrep.R ${yaml}
  while read -r name samyaml; do
    echo -e ".....Running sample $name"
    bash ${ArgenTAG_pipeline}/SCRIPTS/Argentag_Pipeline_V2_main.sh ${samyaml}
  done < ${outdir}/yams.txt
elif [[ ${multi} == "no" ]] ; then
  echo -e "\n....................Preparing single sample pipeline....................\n"
  bash ${ArgenTAG_pipeline}/SCRIPTS/Argentag_Pipeline_V2_main.sh ${yaml}
else 
  echo "ERROR: Multisample paramater must be yes/no"
  exit 1
fi


