#!/bin/bash
####################################################################
#
#A (quite) simple submit script for a one or tow processor job
#
####################################################################
#
# SGE options
#
#Change to the current working directory upon starting of the job
#$ -cwd
#
# Specify the kind of shell script you use, for example, bash
#$ -S /bin/bash
#
# join the error and standard output streams
#$ -j y
#
#
# don't flood myself with e-mail
#$ -m e
#
#
#where the format error go
#$ -e $data_path/log 
#where the format output go
#$ -o $data_path/log
# notify me about pending SIG_STOP and SIG_KILL
#$ -notify
#
### export the default option
source ~/.profile
#name of the job
#$ -N $data_path 
# Specify the array start ,end , step
#$ -t 1-$start:$end:step 
# end of SGE stuff
#########################################################
# now execute my job:
ARRAY=( head W1_R W2_R T1_R T2_R )
###wt1_R1.fastq pp2c_1_R_
#  echo ${ARRAY[$SGE_TASK_ID]}

DIR=/psc/home/zhaocheng/My_project/Dongmeng_SRDX
RE=$DIR/results
DATA=$DIR/data
TMPD=$DIR/tmp_data
CODE=/psc/home/zhaocheng/CLUSTER/code/countsum.pl

/psc/home/zhaocheng/tools/DynamicTrim.pl $DATA/${ARRAY[$SGE_TASK_ID]}1.fastq -h 17 -d $TMP
/psc/home/zhaocheng/tools/DynamicTrim.pl $DATA/${ARRAY[$SGE_TASK_ID]}2.fastq -h 17 -d $TMP
/psc/home/zhaocheng/tools/LengthSort.pl $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.trimmed $TMP/${ARRAY[$SGE_TASK_ID]}2.fastq.trimmed -l 25 -d $TMP

/psc/home/zhaocheng/.local/bin/cutadapt -a AGATCGGAAGAG -f fastq $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.trimmed.paired1 -o $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.cut
/psc/home/zhaocheng/.local/bin/cutadapt -a AGATCGGAAGAG -f fastq $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.trimmed.paired2 -o $TMP/${ARRAY[$SGE_TASK_ID]}2.fastq.cut
/psc/home/zhaocheng/tools/LengthSort.pl $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.cut $TMP/${ARRAY[$SGE_TASK_ID]}2.fastq.cut -l 25 -d $TMP

OUT1="${ARRAY[$SGE_TASK_ID]}1.fastq	`perl $CODE $DATA/${ARRAY[$SGE_TASK_ID]}1.fastq`	`perl $CODE $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.trimmed.paired1`	`perl $CODE $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.cut.paired1`"
OUT2="${ARRAY[$SGE_TASK_ID]}2.fastq	`perl $CODE $DATA/${ARRAY[$SGE_TASK_ID]}2.fastq`	`perl $CODE $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.trimmed.paired2`	`perl $CODE $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.cut.paired2`"
echo "$OUT1
$OUT2" >> $RE/outcome.tmp

mv $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.cut.paired1 $TMPD/${ARRAY[$SGE_TASK_ID]}clean_pair1.fastq
mv $TMP/${ARRAY[$SGE_TASK_ID]}1.fastq.cut.paired2 $TMPD/${ARRAY[$SGE_TASK_ID]}clean_pair2.fastq



# end of job script


