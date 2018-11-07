#!/bin/bash

input="./count0.txt"
p=0
while IFS=' ' read -r f1 f2
do
  if (( $f1 > 1093 )); then
    x=$(( f2 - p ))
    p=$(( p + 1 ))
    sed -i "${x}{d; q}" ./combinfast.txt
    printf 'f2=%s, f1=%s\n' "$f2" "$f1"
  fi
done <"$input"
