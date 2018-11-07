#!/bin/bash

input="./count0.txt"
linelist[1]=0
while IFS=' ' read -r f1 f2
do
  if (( $f1 > 1093 )); then
    linelist[$f2]=1
  fi
done <"$input"
datefile="./combin.txt"
p=1
while read -r line
do
  printf 'processsing line %s ...\n' "$p"
  if [[ ${linelist[$p]} -eq 1 ]]; then
    printf 'delete ...\n' 
  else
    echo "$line" >> "./combinopt.txt"
  fi
  p=$(( p + 1 ))
done <"$datefile"
