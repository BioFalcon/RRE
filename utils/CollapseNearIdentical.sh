#!/bin/bash
#Collapse copies that are near-identical over most of their length (segmental duplications, or the
#same locus sampled twice) into one representative: they are not independent insertions and would
#otherwise dominate the alignment.
#"Near-identical" depends on the repeat: with MaxIdentity=auto, a subsample of the copies is aligned,
#their median pairwise identity measured, and copies more similar than halfway between that and
#identity 1 are collapsed. A number forces a fixed threshold.
#Usage: CollapseNearIdentical.sh <copies.fa> <CPUs> <MaxIdentity|auto>
#The fasta is replaced in place; cd-hit-est clusters are kept as <copies.fa>.clstr

Copies=$1
CPU=$2
MaxIdentity=$3
UtilsDir=$(dirname $(realpath $0))

Before=$(grep -c ">" ${Copies})
[ ${Before} -lt 2 ] && exit 0

if [[ ${MaxIdentity} == "auto" ]];then
    #Median pairwise identity of up to 60 copies
    seqkit shuffle -s 1992 ${Copies} 2> /dev/null | seqkit head -n 60 > ${Copies}.Subsample.fa
    mafft --auto --thread ${CPU} ${Copies}.Subsample.fa > ${Copies}.Subsample.aln.fa 2> /dev/null
    python3 ${UtilsDir}/PairwiseIdentity.py --input ${Copies}.Subsample.aln.fa > ${Copies}.Identity.log
    Threshold=$(awk '$1=="THRESHOLD"{print $2}' ${Copies}.Identity.log)
    Median=$(awk '$1=="MEDIAN_IDENTITY"{print $2}' ${Copies}.Identity.log)
    rm -f ${Copies}.Subsample.fa ${Copies}.Subsample.aln.fa
else
    Threshold=${MaxIdentity}
    Median=NA
fi

#Word size allowed by cd-hit-est for this identity threshold
Word=$(echo ${Threshold} | awk '{if($1>=0.9){print 8}else if($1>=0.88){print 7}else if($1>=0.85){print 6}else{print 5}}')

cd-hit-est \
    -i ${Copies} \
    -o ${Copies}.Collapsed \
    -c ${Threshold} \
    -n ${Word} \
    -G 0 \
    -aS 0.8 \
    -d 0 \
    -M 0 \
    -T ${CPU} > ${Copies}.cdhit.log || exit 0

mv ${Copies}.Collapsed ${Copies}
mv ${Copies}.Collapsed.clstr ${Copies}.clstr
echo -e "\nCopies: median pairwise identity ${Median}; collapsed $(( Before - $(grep -c ">" ${Copies}) )) copies at >= ${Threshold} identity into their representatives"
