#! /bin/bash

mkdir tmp_ep
cp ./01clin_ep_clear.txt ./tmp_ep
cd tmp_ep

awk 'NR%1000==1 { file = "tmp_" sprintf("%05d.txt", NR+999) } { print > file }' 01clin_ep_clear.txt
# split -a 4 -d -l 1000 01clin_ep_clear.txt pval_

rm 01clin_ep_clear.txt

find ./ -name "tmp_*.txt" | parallel '~/Dropbox/LSUHSC/GitSource/bioinfo/gene_match/permtest -f {} -t 916 -p 100000 > {.}.pval'

cat *.pval > ../02clin_ep_permtest.pval

