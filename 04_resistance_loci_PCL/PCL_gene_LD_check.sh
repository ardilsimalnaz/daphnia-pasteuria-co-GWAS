#!/bin/bash

cat > PCL_relevant_genes.tsv <<'GENES'
Pcl21B	317101	317544
Pcl27	317594	318541
Pcl28	318592	319830
Pcl12	1577719	1578984
Pcl13	1579022	1579993
Pcl14	1580131	1581168
Pcl38	1586627	1587892
Pcl37	1587930	1588910
Pcl36	1588930	1589937
GENES

awk 'BEGIN{OFS="\t"}
NR==FNR{
    gene[NR]=$1
    start[NR]=$2
    end[NR]=$3
    n=NR
    next
}
{
    g="between_genes"
    for(i=1;i<=n;i++){
        if($3>=start[i] && $3<=end[i]) g=gene[i]
    }
    print $0,g
}' PCL_relevant_genes.tsv ABC_C_PCL_overlap.tsv \
> ABC_C_PCL_overlap_genes.tsv

awk 'BEGIN{OFS="\t"}
NR==FNR{
    gene[NR]=$1
    start[NR]=$2
    end[NR]=$3
    n=NR
    next
}
{
    g="between_genes"
    for(i=1;i<=n;i++){
        if($3>=start[i] && $3<=end[i]) g=gene[i]
    }
    print $0,g
}' PCL_relevant_genes.tsv D_PCL_overlap.tsv \
> D_PCL_overlap_genes.tsv

plink2 \
  --bfile pasteuria.matched.geno.maf \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 158 \
  --make-bed \
  --out pasteuria.ldids

awk 'NR==FNR {p[$3]=1; next} ($4 in p){print $2}' \
ABC_C_PCL_overlap.tsv pasteuria.ldids.bim \
> ABC_C_PCL.ids

awk 'NR==FNR {p[$3]=1; next} ($4 in p){print $2}' \
D_PCL_overlap.tsv pasteuria.ldids.bim \
> D_PCL.ids

plink2 \
  --bfile pasteuria.ldids \
  --extract ABC_C_PCL.ids \
  --r2-unphased \
  --ld-window-r2 0 \
  --out ABC_C_PCL_LD

plink2 \
  --bfile pasteuria.ldids \
  --extract D_PCL.ids \
  --r2-unphased \
  --ld-window-r2 0 \
  --out D_PCL_LD

echo "----- ABC/C -----"
cut -f5 ABC_C_PCL_overlap_genes.tsv | sort | uniq -c

awk 'NR>1 {
    n++
    if($7 >= 0.8) high++
    if($7 >= 0.5) medium++
    if($7 == 1) perfect++
    sum += $7
}
END {
    print "Total pairs:", n
    print "r2 >= 0.8:", high
    print "r2 >= 0.5:", medium
    print "r2 = 1:", perfect
    print "Mean r2:", sum/n
}' ABC_C_PCL_LD.vcor

echo "----- D -----"
cut -f5 D_PCL_overlap_genes.tsv | sort | uniq -c

awk 'NR>1 {
    n++
    if($7 >= 0.8) high++
    if($7 >= 0.5) medium++
    if($7 == 1) perfect++
    sum += $7
}
END {
    print "Total pairs:", n
    print "r2 >= 0.8:", high
    print "r2 >= 0.5:", medium
    print "r2 = 1:", perfect
    print "Mean r2:", sum/n
}' D_PCL_LD.vcor
