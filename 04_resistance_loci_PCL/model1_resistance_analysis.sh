#!/bin/bash

awk -F'\t' 'BEGIN{OFS="\t"}

NR==1 {
    print $0,"DAPH_POS","RESISTANCE_LOCUS"
    next
}

{
    x=$3
    sub(/^[^:]*:/,"",x)
    sub(/[A-Za-z].*/,"",x)
    pos=x+0

    locus=""

    if ($2=="000011F|quiver" && pos>=2170278 && pos<=2357097)
        locus="ABC_C"

    else if ($2=="000018F|quiver" && pos>=1801000 && pos<=1804000)
        locus="D"

    else if ($2=="000067F|quiver" && pos>=233835 && pos<=388722)
        locus="E"

    if (locus!="")
        print $0,pos,locus
}
' model1_KEEP/model1_bonferroni_significant.tsv \
> model1_resistance_loci_overlap.tsv


awk 'END{print "Total overlapping associations:", NR-1}' \
model1_resistance_loci_overlap.tsv


awk -F'\t' '
NR>1 {count[$NF]++}
END {
    for (locus in count)
        print locus, count[locus]
}
' model1_resistance_loci_overlap.tsv


for locus in ABC_C D E
do
    echo -n "$locus unique Daphnia variants: "
    awk -F'\t' -v locus="$locus" \
    'NR>1 && $NF==locus {print $3}' \
    model1_resistance_loci_overlap.tsv | sort -u | wc -l
done


for locus in ABC_C D E
do
    echo -n "$locus unique Pasteuria phenotypes: "
    awk -F'\t' -v locus="$locus" \
    'NR>1 && $NF==locus {print $1}' \
    model1_resistance_loci_overlap.tsv | sort -u | wc -l
done


awk -F'\t' 'NR>1 && $NF=="ABC_C" {print $1}' \
model1_resistance_loci_overlap.tsv | sort -u \
> ABC_C_PAST_ids.txt

awk -F'\t' 'NR>1 && $NF=="D" {print $1}' \
model1_resistance_loci_overlap.tsv | sort -u \
> D_PAST_ids.txt


awk 'NR==FNR {keep[$1]=1; next} ($1 in keep)' \
ABC_C_PAST_ids.txt pasteuria_pheno_map.tsv \
> ABC_C_pasteuria_variants.tsv

awk 'NR==FNR {keep[$1]=1; next} ($1 in keep)' \
D_PAST_ids.txt pasteuria_pheno_map.tsv \
> D_pasteuria_variants.tsv


awk -F'\t' 'BEGIN{OFS="\t"}
{
    x=$2
    sub(/^[^:]*:/,"",x)
    sub(/[A-Za-z].*/,"",x)
    print $1,$2,x
}' ABC_C_pasteuria_variants.tsv \
> ABC_C_pasteuria_coordinates.tsv

awk -F'\t' 'BEGIN{OFS="\t"}
{
    x=$2
    sub(/^[^:]*:/,"",x)
    sub(/[A-Za-z].*/,"",x)
    print $1,$2,x
}' D_pasteuria_variants.tsv \
> D_pasteuria_coordinates.tsv


cat > PCL_triplets.tsv <<'EOF'
PCL_triplet_1_Pcl6_7_8	7278	10679
PCL_triplet_2_Pcl23_22_21A	86443	90028
PCL_triplet_3_Pcl21B_27_28	317101	319830
PCL_triplet_4_Pcl15_16_17	1199952	1203511
PCL_triplet_5_Pcl24_25_26	1517919	1521372
PCL_triplet_6_Pcl9_10_11	1571213	1574559
PCL_triplet_7_Pcl12_13_14	1577719	1581168
PCL_triplet_8_Pcl38_37_36	1586627	1589937
EOF


awk -F'\t' 'BEGIN{OFS="\t"}
NR==FNR {
    name[++n]=$1
    start[n]=$2
    end[n]=$3
    next
}
{
    pos=$3
    for(i=1;i<=n;i++)
        if(pos>=start[i] && pos<=end[i])
            print $0,name[i]
}' PCL_triplets.tsv ABC_C_pasteuria_coordinates.tsv \
> ABC_C_PCL_overlap.tsv


awk -F'\t' 'BEGIN{OFS="\t"}
NR==FNR {
    name[++n]=$1
    start[n]=$2
    end[n]=$3
    next
}
{
    pos=$3
    for(i=1;i<=n;i++)
        if(pos>=start[i] && pos<=end[i])
            print $0,name[i]
}' PCL_triplets.tsv D_pasteuria_coordinates.tsv \
> D_PCL_overlap.tsv


wc -l ABC_C_PCL_overlap.tsv D_PCL_overlap.tsv


for file in ABC_C_PCL_overlap.tsv D_PCL_overlap.tsv
do
    echo
    echo "===== $file ====="

    awk -F'\t' '{count[$4]++}
    END {
        for(p in count)
            print p, count[p]
    }' "$file" | sort -k2,2nr
done


awk -F'\t' 'BEGIN{OFS="\t"}
NR==FNR {
    name[++n]=$1
    start[n]=$2
    end[n]=$3
    next
}
{
    pos=$3
    min=-1
    nearest=""

    for(i=1;i<=n;i++) {
        if(pos>=start[i] && pos<=end[i])
            d=0
        else if(pos<start[i])
            d=start[i]-pos
        else
            d=pos-end[i]

        if(min==-1 || d<min) {
            min=d
            nearest=name[i]
        }
    }

    print $0,nearest,min
}' PCL_triplets.tsv ABC_C_pasteuria_coordinates.tsv \
> ABC_C_nearest_PCL.tsv


awk -F'\t' 'BEGIN{OFS="\t"}
NR==FNR {
    name[++n]=$1
    start[n]=$2
    end[n]=$3
    next
}
{
    pos=$3
    min=-1
    nearest=""

    for(i=1;i<=n;i++) {
        if(pos>=start[i] && pos<=end[i])
            d=0
        else if(pos<start[i])
            d=start[i]-pos
        else
            d=pos-end[i]

        if(min==-1 || d<min) {
            min=d
            nearest=name[i]
        }
    }

    print $0,nearest,min
}' PCL_triplets.tsv D_pasteuria_coordinates.tsv \
> D_nearest_PCL.tsv


for file in ABC_C_nearest_PCL.tsv D_nearest_PCL.tsv
do
    echo
    echo "===== $file ====="

    awk -F'\t' '
    {
        d=$NF

        if(d==0)
            inside++
        else if(d<=1000)
            kb1++
        else if(d<=5000)
            kb5++
        else if(d<=10000)
            kb10++
        else
            far++
    }
    END {
        print "Inside PCL triplet:", inside+0
        print "Within 1 kb:", kb1+0
        print "1-5 kb:", kb5+0
        print "5-10 kb:", kb10+0
        print ">10 kb:", far+0
    }' "$file"
done


for file in ABC_C_nearest_PCL.tsv D_nearest_PCL.tsv
do
    echo
    echo "===== $file ====="

    awk -F'\t' '{count[$4]++}
    END {
        for(p in count)
            print p, count[p]
    }' "$file" | sort -k2,2nr
done