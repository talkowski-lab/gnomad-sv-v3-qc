#!/bin/bash

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Collects QC data for SV VCF output by SV pipeline
# Cohort-level benchmarking comparisons to external dataset

set -e

###USAGE
usage(){
cat <<EOF

usage: collectQC.external_benchmarking.sh [-h] [-r REF] STATS SVTYPES CONTIGS BENCHDIR OUTDIR

Helper tool to collect cohort-level bechmarking data for a VCF output by sv-pipeline vs. an external dataset

Positional arguments:
  STATS           VCF stats file generated by collectQC.vcf_wide.sh
  SVTYPES         List of SV types to evaluate. Two-column, tab-delimited file.
                  First column: sv type. Second column: HEX color for sv type.
  CONTIGS         List of contigs to evaluate
  BENCHDIR        Directory containing benchmark archives
  OUTDIR          Output directory for all QC data

Optional arguments:
  -h  HELP        Show this help message and exit
  -q  QUIET       Silence all status updates

EOF
}


###PARSE ARGS
QUIET=0
while getopts ":qh" opt; do
	case "$opt" in
		h)
			usage
			exit 0
			;;
    q)
      QUIET=1
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
STATS=$1
SVTYPES=$2
CONTIGS=$3
BENCHDIR=$4
OUTDIR=$5


###PROCESS ARGS & OPTIONS
#Check for required input
if [ -z ${STATS} ]; then
  echo -e "\nERROR: input VCF stats file not specified\n"
  usage
  exit 0
fi
if ! [ -s ${STATS} ]; then
  echo -e "\nERROR: input VCF stats file either empty or not found\n"
  usage
  exit 0
fi
#Check for required list of contigs
if [ -z ${CONTIGS} ]; then
  echo -e "\ncollectQC.external_benchmarking.sh ERROR: input contig list not specified\n"
  usage
  exit 1
elif ! [ -s ${CONTIGS} ]; then
  echo -e "\ncollectQC.external_benchmarking.sh ERROR: input contig list either empty or not found\n"
  usage
  exit 1
fi
if [ -z ${OUTDIR} ]; then
  echo -e "\nERROR: output directory not specified\n"
  usage
  exit 0
fi


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###PREP INPUT FILES
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Preparing input files for external benchmarking"
fi
#Prep directories
QCTMP=`mktemp -d`
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
if ! [ -e ${OUTDIR}/data ]; then
  mkdir ${OUTDIR}/data
fi
#Gather SV types to process
cut -f1 ${SVTYPES} | sort | uniq > ${QCTMP}/svtypes.txt
# Gather list of external BEDs to compare
for bed in $( find ${BENCHDIR} -name "*.bed.gz" ); do
  prefix=$( basename $bed | sed 's/\.bed\.gz//g' )
  echo -e "${prefix}\t${bed}"
done > comparator_beds.tsv


###GATHER EXTERNAL BENCHMARKING
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Starting external benchmarking"
fi
while read prefix bed; do
  #Print status
  if [ ${QUIET} == 0 ]; then
    echo -e "$( date ) - VCF QC STATUS: Benchmarking samples in ${prefix}"
  fi
  ${BIN}/compare_callsets.sh \
    -O ${QCTMP}/${prefix}.overlaps.bed \
    -p ${prefix}_Benchmarking_SV \
    ${STATS} \
    ${bed} \
    ${CONTIGS}
  cp ${QCTMP}/${prefix}.overlaps.bed ${OUTDIR}/data/
  bgzip -f ${OUTDIR}/data/${prefix}.overlaps.bed
  tabix -f ${OUTDIR}/data/${prefix}.overlaps.bed.gz
done < comparator_beds.tsv


###CLEAN UP
rm -rf ${QCTMP}

