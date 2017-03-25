#!/bin/bash
# Reading in fastq.list
while read line; do
  file_name=$(basename $line);
  base_name="${file_name%.*}";
  echo -en $base_name'\t';
  line_count=`wc -l $line | awk {'print $1'}`; echo "$line_count/4" | bc;
done < $1
