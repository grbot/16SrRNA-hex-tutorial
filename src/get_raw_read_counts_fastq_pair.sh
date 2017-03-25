#!/bin/bash
# Reading in sid.fastq_pair.list
while read line; do 
  sid=`echo $line | awk '{print $1}'`; echo -en $sid'\t'; 
  read1=`echo $line | awk '{print $2}'`; line_count=`wc -l $read1 | awk {'print $1'}`;read1_count=`echo "$line_count/4" | bc`; echo -en $read1_count'\t';
  read2=`echo $line | awk '{print $3}'`; line_count=`wc -l $read2 | awk {'print $1'}`;read2_count=`echo "$line_count/4" | bc`; echo -e $read2_count;
done < $1


