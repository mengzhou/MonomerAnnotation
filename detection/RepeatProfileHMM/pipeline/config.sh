#!/bin/bash
# Set number of threads for parallelization
NTHREAD=1
# Set path to dependencies here
# 1. path to Tandem Repeats Finder executable binary
TRF=
# 2. path to MUSCLE executable binary
MUSCLE=
# 3. path to HMMER *directory* of binaries. E.g.: /home/user/HMMER/binaries
HMMER=
# 4. path to bedtools *directory* of binaries
BEDTOOLS=
# 5. path to RepeatProfileHMM *root* directory
RPHMM=
# 6. set path to reference genome, must be a single FASTA file
REF=

# validation of paths. no need to change things below
REV_COMP=$RPHMM/utils/revcomp_fa.py
GET_MEDIAN=$RPHMM/utils/get_median.awk
NUMBERING=$RPHMM/utils/numbering.awk
for i in $TRF $MUSCLE $HMMER $NHMMER $BEDTOOLS $RPHMM $REF $REV_COMP
do
  if [ ! -e $i ]
  then
    echo "Path $i does not exist! Please check config.sh"
    exit
  fi
done
