#!/bin/bash

##-------------------------------------------------------------------------------------------------
XLSFILE=BRCA.merged_only_clinical_clin_format
SRCFILE=BRCA.merged_only_clinical_clin_format

libreoffice --headless --convert-to csv ${XLSFILE}.xlsx --outdir ./

## Method 1
# cut -d ',' -f1 PRAD_clinical.csv > temp

## Methos 1.1
# awk -F' ' 'NR==FNR{a[$1]; next}; $1 in a' temp BRCA.merged_only_clinical_clin_format.txt > BRCA_clinical.txt

## Method 1.2
# grep -w -f temp BRCA.merged_only_clinical_clin_format.txt > BRCA_clinical.txt

# rm temp



## Method 2
cut -d ',' -f1 ${XLSFILE}.csv | grep -w -f - ${SRCFILE}.txt > ${SRCFILE}_cut.txt

rm ${XLSFILE}.csv


##-------------------------------------------------------------------------------------
CLIN_EXP=clin_exp_1.txt

mkdir ./tmp

tail -n+31 00clin_ep.txt > 01clin_ep.txt

join -1 1 -2 1 -t $'\t' <(sort -k 1,1 ../../RefGene/Gene_location.txt) <(paste <(cut -f1 01clin_trip.txt | awk -F '\t' 'BEGIN{OFS="\t";} { sub(/\..+$/, "", $1) }1') 01clin_trip.txt | sort -k 1,1) | cut -f5- > 01clin_trip_syb.txt

head -29 00clin_ep.txt | tail -1 > 01clin_ep_index.txt
split -l 21728 01clin_ep.txt
cat <(tail -n+15 02clin_ep_xaa_multest_parallel.res) <(tail -n+15 02clin_ep_xab_multest_parallel.res) | nl | cut -f 1,3- > 02clin_ep_multest_parallel.txt
join -1 2 -2 1 <(sort -k 2,2 02clin_ep_multest_parallel.txt) <(sort -k 1,1 01clin_ep.txt) -t $'\t' | sort -g -k 2.2 > 03clin_ep_p-val.txt
awk -F '\t' 'BEGIN{OFS="\t";} { sub(/\..+$/, "", $1) }1' 03clin_ep_p-val.txt > 04clin_ep_p-val.txt
join -1 1 -2 1 -t $'\t' <(sort -k 1,1 ../../RefGene/Gene_location.txt) <(sort -k 1,1 04clin_ep_p-val.txt) | cut -f 1,5- | sort -g -k 4.4 > 05clin_ep_p-val.txt
head -200 05clin_ep_p-val.txt | cut -f 2,8- > 06clin_ep_p-val.txt 

join -1 2 -2 1 <(sort -k 2,2 02clin_trip_multest_parallel.txt) <(sort -k 1,1 01clin_trip.txt) -t $'\t' | sort -g -k 2,2 > 03clin_trip_p-val.txt
awk -F '\t' 'BEGIN{OFS="\t";} { sub(/\..+$/, "", $1) }1' 03clin_trip_p-val.txt > 04clin_trip_p-val.txt
join -1 1 -2 1 -t $'\t' <(sort -k 1,1 ../../RefGene/Gene_location.txt) <(sort -k 1,1 04clin_trip_p-val.txt) | cut -f 1,5- | sort -g -k 4,4 > 05clin_trip_p-val.txt,
head -200 05clin_trip_p-val.txt | cut -f 2,8- > 06clin_trip_p-val.txt 
cat <(awk 'NR==1 || NR==9 || NR==30' 00clin_trip.txt) <(cut -f 2,8- 05clin_trip_p-val.txt) > 06clin_trip_p-val_whole.txt

##----------------------------------------------------------------------------
#for permutation test
cut -f 2-83,85-491,493-504,506-560,562-580,582- 01clin_ep.txt > 01clin_ep_clear.txt
~/Dropbox/LSUHSC/GitSource/bioinfo/gene_match/permtest -f 01clin_ep_clear.txt -t 916 -p 100000 > 02clin_ep_permtest.pval &
# or...
./perm_test_ep.sh&
join -1 1 -2 1 -t $'\t' <(sort -k 1,1 ../../RefGene/Gene_location.txt) <(paste <(cut -f1 01clin_ep.txt | awk -F '\t' 'BEGIN{OFS="\t";} { sub(/\..+$/, "", $1) }1') 02clin_ep_permtest.pval | sort -k 1,1) | cut -f 1,5- | sort -k3,3gr | awk '{print $1,$2,$3,-10*log($3)/log(10)}' OFS=$'\t' > 05clin_ep_permtest_pval.txt

cut -f 2-49,51- 01clin_trip.txt > 01clin_trip_clear.txt
~/Dropbox/LSUHSC/GitSource/bioinfo/gene_match/permtest -f 01clin_trip_clear.txt -t 110 -p 100000 > 02clin_trip_permtest.pval &
# or...
./perm_test_trip.sh&
# paste <(cut -f1 01clin_trip.txt | awk -F '\t' 'BEGIN{OFS="\t";} { sub(/\..+$/, "", $1) }1') 02clin_trip_permtest.pval | sort -k 1,1  > 04clin_trip_permtest_pval.txt
join -1 1 -2 1 -t $'\t' <(sort -k 1,1 ../../RefGene/Gene_location.txt) <(paste <(cut -f1 01clin_trip.txt | awk -F '\t' 'BEGIN{OFS="\t";} { sub(/\..+$/, "", $1) }1') 02clin_trip_permtest.pval | sort -k 1,1) | cut -f 1,5- | sort -k3,3gr | awk '{print $1,$2,$3,-10*log($3)/log(10)}' OFS=$'\t' > 05clin_trip_permtest_pval.txt


##-----------------------------------------------------------------------------------------------------------------
#field 1: Hugo_Symbol; 5: Chromosome; 9: Variant_Classification; 10: Variant_Type; 40: cDNA_Change
mkdir result
find ./ -name "TCGA*.txt" -printf "%f\n" | parallel 'VAR={.}; VAR=${VAR/.*/}; tail -n+5 {} | cut -f1,10 > ./result/${VAR}.txt'
cd result
find ./ -name "TCGA*.txt" -printf "%f\n" | parallel "sed -i '1s/^/Hugo_Symbol\t{.}\n/' {}"
~/Dropbox/LSUHSC/GitSource/MyScripts/full-outer-join.awk HEADER=1 TCGA*.txt | sed 's/NULL//g' > Variant_Type.txt
rm TCGA*.txt

grep -o -n -E 'SNP|DEL|INS' Variant_Type.txt | cut -d : -f 1 | uniq -c | awk 'NR==FNR{a[$2]=$1; next}; {print a[FNR]"\t"$0}' - Variant_Type.txt | sed '1s/^/ALL#/' | awk 'NR<2{print $0;next}{print $0 | "sort -k1rn"}' > Variant_Type_ALL.txt
grep -o -n 'SNP' Variant_Type_ALL.txt | cut -d : -f 1 | uniq -c | awk 'NR==FNR{a[$2]=$1; next}; {print a[FNR]"\t"$0}' - Variant_Type_ALL.txt | sed '1s/^/SNP#/' > Variant_Type_SNP.txt
grep -o -n -E 'DEL|INS' Variant_Type_ALL.txt | cut -d : -f 1 | uniq -c | awk 'NR==FNR{a[$2]=$1; next}; {print a[FNR]"\t"$0}' - Variant_Type_SNP.txt | sed '1s/^/INDEL#/' > Variant_Type_SNP_INDEL.txt
grep -o -n -E 'DEL|INS' Variant_Type.txt | cut -d : -f 1 | uniq -c | awk 'NR==FNR{a[$2]=$1; next}; {print a[FNR]"\t"$0}' - Variant_Type.txt | sed '1s/^/INDEL#/' | awk 'NR<2{print $0;next}{print $0 | "sort -k1rn"}' > Variant_Type_INDEL.txt
# or...
awk -F $'\t' 'NR<2{print $0;next}{s=0;d=0;n=0;o=0;for(i=2;i<=NF;i++)if($i~/SNP/)s++;else if($i~/DEL/)d++;else if($i~/INS/)n++;else if($i!="")o++; print s,d,n,o,s+d+n+o,100*(s+d+n+o)/(NF-1),$0 | "sort -k5rn"}' OFS=$'\t' Variant_Type.txt | sed '1s/^/SNP\tDEL\tINS\tOTHER\tALL\tfreq(%)\t/' > Variant_Type_Count.txt

~/Dropbox/LSUHSC/GitSource/MyScripts/full-outer-join.awk <(awk 'NR>1{print $7,$5,$6}' OFS=$'\t' Variant_Type_Count.txt) <(cut -f2,4 ../../../pomelo_tripnega_TPM/05clin_trip_p-val.txt) | sort -k2rn | awk '{print $1,$2,$3,-10*log($4)/log(10),$4}' OFS=$'\t' | sed -e 's/NULL//' -e '1iHugo_sym\tMut_count\tMut_freq(%)\tp_Score\tp-value' > mut_exp.txt
# add... permutation test
~/Dropbox/LSUHSC/GitSource/MyScripts/full-outer-join.awk mut_exp.txt <(cut -f2,4 ../../../pomelo_tripnega_TPM/05clin_trip_permtest_pval.txt) | sort -k2rn | sed -e 's/NULL//' -e '1iHugo_sym\tMut_count\tMut_freq(%)\t-10log(Ttest)\tTtest_pval\t-10log(Ptest)' > mut_exp_permtest.txt

# for ep
~/Dropbox/LSUHSC/GitSource/MyScripts/full-outer-join.awk mut_exp.txt <(cut -f2,4 ../../../pomelo_tripnega_TPM/05clin_ep_permtest_pval.txt) | sort -k2rn | sed -e 's/NULL//' -e '1iHugo_sym\tMut_count\tMut_freq(%)\t-10log(Ttest)\tTtest_pval\t-10log(Ptest)' > mut_exp_permtest.txt

##--------------------------------------------------------------------------------
# awk '($2>0 && $3>30 && $3!="inf")' mut_exp.txt > mut_exp_filtered.txt
awk -F $'\t' '$1=$1' gwas.txt | sort | uniq | awk 'BEGIN{print "Hugo_sym\tGWAS";} {print $1,1}' OFS=$'\t' > gwas_uniq.txt
~/Dropbox/LSUHSC/GitSource/MyScripts/full-outer-join.awk HEADER=1 ./trip/result/mut_exp.txt ./nontrip/result/mut_exp.txt ./gwas_uniq.txt | cut -f1,2,3,4,6,7,8,10 | sed -e '1s/^.*$/Hugo_sym\tTrip_Mu\tTrip_Mu(%)\tTrip_pscore\tEp_Mu\tEp_Mu(%)\tEp_pscore\tGWAS/' -e 's/NULL//g' | awk -F $'\t' '($3>0 || $4>30 || $6>0 || $7>30 || $8>0)' | cut -f1,3,4,6,7,8 > mut_exp.txt
# or ... for permutation test
~/Dropbox/LSUHSC/GitSource/MyScripts/full-outer-join.awk HEADER=1 ./trip/result/mut_exp_permtest.txt ./nontrip/result/mut_exp_permtest.txt ./gwas_uniq.txt | cut -f1,3,6,8,11,12 | sed -e 's/NULL//g' | awk -F $'\t' '($2>0 || $3>33 || $4>0 || $5>33 || $6>0)' | sed -e '1s/^.*$/Hugo_sym\tTrip_Mu(%)\ttrip-10log(Ptest)\tEp_Mu(%)\tEp-10log(Ptest)\tGWAS/' > mut_exp_permtest.txt


awk -F $'\t' '(NR>1 && $2>0){print $1}' mut_exp.txt > trip_mu0.txt
awk -F $'\t' '(NR>1 && $4>0){print $1}' mut_exp.txt > ep_mu0.txt

awk -F $'\t' '(NR>1 && $3>30){print $1}' mut_exp.txt > trip_exp30.txt
awk -F $'\t' '(NR>1 && $5>30){print $1}' mut_exp.txt > ep_exp30.txt

awk -F $'\t' '(NR>1 && ($2>0 || $4>0)){print $1,$2,$4}' OFS=$'\t' mut_exp.txt > mut_freq.txt

##---------------------------
#permutation test
awk -F $'\t' '(NR>1 && $2>0){print $1}' mut_exp_permtest.txt > trip_mu0_permtest.txt
awk -F $'\t' '(NR>1 && $4>0){print $1}' mut_exp_permtest.txt > ep_mu0_permtest.txt

awk -F $'\t' '(NR>1 && $3>33){print $1}' mut_exp_permtest.txt > trip_exp33_permtest.txt
awk -F $'\t' '(NR>1 && $5>33){print $1}' mut_exp_permtest.txt > ep_exp33_permtest.txt


join --header -1 1 -2 1 -t $'\t' ../../RefGene/cytoband.txt mut_exp_permtest.txt | cut -f1-4,7 | awk -F $'\t' '(NR<2){next}{print $1,$2,10^(-$4/10),int(0.5+$3*99/100),$5,$4}' OFS=$'\t' | sed -e '1iGenes\tPosition\tp-value\tSomatic\tGWAS' > 33.txt
awk -F $'\t' '(NR>1 && $4>0){print $1}' 33.txt > trip_mu0_permtest-33.txt
awk -F $'\t' '(NR>1 && $6>33){print $1}' 33.txt > trip_exp33_permtest-33.txt
awk -F $'\t' '(NR>1 && $5==1){print $1}' 33.txt > trip_gwas-33.txt
cut -f1-4 33.txt | awk -F $'\t' '(NR>1 && $4>0){print}' | sort -t $'\t' -k4,4rn | sed -e '1iGenes\tPosition\tp-value\tSomatic' > trip_Somatic.txt
awk -F $'\t' '(NR>1 && $5==1){print}' 33.txt | sort -t $'\t' -k4,4rn | sed -e '1iGenes\tPosition\tp-value\tSomatic\tGWAS' > trip_Somatic_GWAS.txt

join --header -1 1 -2 1 -t $'\t' ../../RefGene/cytoband.txt mut_exp_permtest.txt | cut -f1,2,5,6,7 | awk -F $'\t' '(NR<2){next}{print $1,$2,10^(-$4/10),int(0.5+$3*7.68),$5}' OFS=$'\t' | sed -e '1iGenes\tPosition\tp-value\tSomatic\tGWAS' > 44.txt
cut -f1-4 44.txt | awk -F $'\t' '(NR>1 && $4>0){print}' | sort -t $'\t' -k4,4rn | sed -e '1iGenes\tPosition\tp-value\tSomatic' > ep_Somatic.txt
awk -F $'\t' '(NR>1 && $5==1){print}' 44.txt | sort -t $'\t' -k4,4rn | sed -e '1iGenes\tPosition\tp-value\tSomatic\tGWAS' > ep_Somatic_GWAS.txt

join --header -1 1 -2 1 -t $'\t' ../../RefGene/cytoband.txt <(awk -F $'\t' '(NR>1 && $3>33 && $5>33){print}' mut_exp_permtest.txt) | cut -f1-7 | awk -F $'\t' '(NR<2){next}{print $1,$2,10^(-$4/10),10^(-$6/10),int(0.5+$3*99/100),int(0.5+$5*7.68),$7}' OFS=$'\t' | sort -t $'\t' -k7,7rn -k6,6rn -k5,5rn | sed -e '1iGenes\tPosition\ttrip_pvalue\tEP_pvalue\ttrip_Somatic\tEp_Somatic\tGWAS' > trip_ep_Somatic_GWAS.txt


## 8/22/2018
# vim: Insert line numbers
# The :s command can be used to insert line numbers before each line:
# The pattern ^ matches the start of every line, 
# and the replacement \= is the result of evaluating the following expression. 
# That expression uses printf() to format the number of the current line
:%s/^/\=printf("%d\t", line('.'))


╭─jiande@BIG-L08 ~/prjs/BreastTumor/Tumor-txt  
╰─$ grep -wFf TNBClist.txt <(grep -E 'file_name|case_id' ../files.2017-12-12T19_25_53.669792.json | awk -F':' '{print $2}' | sed -e 's/"//g;s/,//' | sed -e 'N;s/\n/\t/') > TNBClist_patientid.txt
