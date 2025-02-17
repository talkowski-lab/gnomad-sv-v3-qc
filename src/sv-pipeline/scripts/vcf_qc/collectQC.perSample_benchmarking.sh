#!/bin/bash

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Collects QC data for SV VCF output by SV pipeline
# Per-sample benchmarking comparisons to external dataset

set -e

###USAGE
usage(){
cat <<EOF

usage: collectQC.perSample_benchmarking.sh [options] VCFSTATS SAMPLES CONTIGS PERSAMPDIR.tar.gz SET2.tar.gz OUTDIR

Helper tool to collect per-sample bechmarking data for a VCF output by sv-pipeline vs. an external dataset

Positional arguments:
  VCFSTATS               VCF stats file generated by collectQC.vcf_wide.sh
  SAMPLES                List of samples to be considered
  CONTIGS                List of contigs to evaluate
  PERSAMPDIR.tar.gz      Archive of directory with per-sample variant IDs
  SET2.tar.gz            Ground truth SV callset to use for benchmarking
  OUTDIR                 Output directory for per-sample comparison results

Optional arguments:
  -h  HELP        Show this help message and exit
  -d  DISTANCE    Maximum distance between breakpoints during comparisons (default: 250bp)
  -p  PREFIX      Prefix for benchmarking variant IDs (default: Benchmarking_SV)
  -q  QUIET       Silence all status updates

Notes:
  1) SET2 tarball is expected to contain one BED3+ formatted file per sample
  2) SET2 files must have exactly six columns as follows, in order:
     chr, start, end, SV type, SV size, allele frequency
  3) Assumes all per-sample files SET2 follow the same file naming convention.
     Filename scheme: [sample_ID].[callset_name].[other_suffixes].SV_calls.bed.gz

EOF
}


###PARSE ARGS
PREFIX="Benchmarking_SV"
DIST=250
QUIET=0
while getopts ":r:p:d:qh" opt; do
	case "$opt" in
		h)
			usage
			exit 1
			;;
    p)
      PREFIX=${OPTARG}
      ;;
    d)
      DIST=${OPTARG}
      ;;
    q)
      QUIET=1
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
VCFSTATS=$1
SAMPLES=$2
CONTIGS=$3
PERSAMPDIR_TAR=$4
SET2=$5
OUTDIR=$6


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###PROCESS ARGS & OPTIONS
#Check for required input VCFSTATS
if [ -z ${VCFSTATS} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input variant stats file not specified\n"
  usage
  exit 1
elif ! [ -s ${VCFSTATS} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input variant stats file either empty or not found\n"
  usage
  exit 1
fi
#Check for required list of samples
if [ -z ${SAMPLES} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input sample list not specified\n"
  usage
  exit 1
elif ! [ -s ${SAMPLES} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input sample list either empty or not found\n"
  usage
  exit 1
fi
#Check for required list of contigs
if [ -z ${CONTIGS} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input contig list not specified\n"
  usage
  exit 1
elif ! [ -s ${CONTIGS} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input contig list either empty or not found\n"
  usage
  exit 1
fi
#Check for required input PERSAMPDIR tarball
if [ -z ${PERSAMPDIR_TAR} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input per-sample variant lists tarball not specified\n"
  usage
  exit 1
elif ! [ -s ${PERSAMPDIR_TAR} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input per-sample variant lists tarball either empty or not found\n"
  usage
  exit 1
fi
#Check for required input SET2
if [ -z ${SET2} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input SET2.tar.gz not specified\n"
  usage
  exit 1
elif ! [ -s ${SET2} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: input SET2.tar.gz either empty or not found\n"
  usage
  exit 1
fi
#Checks output directory
if [ -z ${OUTDIR} ]; then
  echo -e "\ncollectQC.perSample_benchmarking.sh ERROR: output directory not specified\n"
  usage
  exit 1
fi
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
#Checks prefix
if [ -z ${PREFIX} ]; then
  PREFIX="Benchmarking_SV"
fi

###PREP INPUT FILES
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - collectQC.perSample_benchmarking.sh STATUS: Preparing input files..."
fi
#Make temporary directory
OVRTMP=`mktemp -d`
#Untar per-sample directory
mkdir ${OVRTMP}/persamp/
tar -xzvf ${PERSAMPDIR_TAR} \
  --directory ${OVRTMP}/persamp/
#Create matching SET1.tar.gz file from VCFSTATS
mkdir ${OVRTMP}/SET1_calls/
while read ID; do
  VIDlist=$( find ${OVRTMP}/persamp/ -name "${ID}.VIDs_genotypes.txt.gz" )
  if [ ! -z ${VIDlist} ] && [ -s ${VIDlist} ]; then
    echo $ID
    zcat ${VCFSTATS} | head -n1 > \
    ${OVRTMP}/SET1_calls/${ID}.SET1.SV_calls.bed
    zcat ${VIDlist} | cut -f1 | fgrep -wf - <( zcat ${VCFSTATS} ) >> \
    ${OVRTMP}/SET1_calls/${ID}.SET1.SV_calls.bed
    bgzip -f ${OVRTMP}/SET1_calls/${ID}.SET1.SV_calls.bed
    tabix -f ${OVRTMP}/SET1_calls/${ID}.SET1.SV_calls.bed.gz
  fi
done < ${SAMPLES}
tar -czvf ${OVRTMP}/SET1_calls.tar.gz ${OVRTMP}/SET1_calls


###RUN COMPARISONS
${BIN}/compare_callsets_perSample.sh \
  -d ${DIST} \
  -p ${PREFIX} \
  ${OVRTMP}/SET1_calls.tar.gz \
  ${SET2} \
  ${SAMPLES} \
  ${CONTIGS} \
  ${OUTDIR}


###CLEAN UP
rm -rf ${OVRTMP}

