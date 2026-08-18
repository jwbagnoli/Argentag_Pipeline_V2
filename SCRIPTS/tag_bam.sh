#!/bin/bash

input=$1
output=$2
bc_tag=$3
umi_tag=$4
threads=$5

samtools view -@ ${threads}  -h ${input} | \
mawk -v bc_tag=${bc_tag} -v umi_tag=${umi_tag} '
BEGIN { OFS="\t" }
{
    if ($0 ~ /^@/) {
        print $0
    } else {
        qname = $1
        split(qname, a, "#")
        split(a[1], b, "_")

        barcode = b[1]
        umi = b[2]

        if (barcode != "" && umi != "") {
            print $0, bc_tag ":Z:" barcode, umi_tag ":Z:" umi
        } else {
            print $0
        }
    }
}' | \
samtools view -@ ${threads} -b -o ${output}




