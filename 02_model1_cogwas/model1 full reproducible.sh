#!/usr/bin/env bash
set -euo pipefail
cd ~/coGWAS_paperstyle/model1_from_onevariant
pwd
/home/go69lef/bin/plink2 --version
ls -lh daphnia_input.vcf pasteuria_input.vcf.gz
file daphnia_input.vcf pasteuria_input.vcf.gz
/home/go69lef/bin/plink2 --vcf daphnia_input.vcf --double-id --allow-extra-chr --max-alleles 2 --make-bed --out daphnia
wc -l daphnia.bim daphnia.fam
/home/go69lef/bin/plink2 --vcf pasteuria_input.vcf.gz --double-id --allow-extra-chr --max-alleles 2 --make-bed --out pasteuria
wc -l pasteuria.bim pasteuria.fam
/home/go69lef/bin/plink2 --bfile daphnia --allow-extra-chr --geno 0.10 --make-bed --out daphnia.geno
wc -l daphnia.geno.bim daphnia.geno.fam
/home/go69lef/bin/plink2 --bfile daphnia.geno --allow-extra-chr --mind 0.10 --make-bed --out daphnia.geno.mind
wc -l daphnia.geno.mind.bim daphnia.geno.mind.fam
/home/go69lef/bin/plink2 --bfile daphnia.geno.mind --allow-extra-chr --maf 0.05 --make-bed --out daphnia.geno.mind.maf
wc -l daphnia.geno.mind.maf.bim daphnia.geno.mind.maf.fam
/home/go69lef/bin/plink2 --bfile daphnia.geno.mind.maf --allow-extra-chr --king-cutoff 0.354 --make-bed --out daphnia.king354
wc -l daphnia.king354.bim daphnia.king354.fam daphnia.king354.king.cutoff.out.id
/home/go69lef/bin/plink2 --bfile daphnia.king354 --allow-extra-chr --set-all-var-ids '@:#$r,$a' --new-id-max-allele-len 243 --indep-pairwise 200kb 1 0.2 --out daphnia.prune
wc -l daphnia.prune.prune.in daphnia.prune.prune.out
/home/go69lef/bin/plink2 --bfile daphnia.king354 --allow-extra-chr --set-all-var-ids '@:#$r,$a' --new-id-max-allele-len 243 --extract daphnia.prune.prune.in --pca 10 --out daphnia.pca
wc -l daphnia.pca.eigenvec
/home/go69lef/bin/plink2 --bfile daphnia.king354 --allow-extra-chr --set-all-var-ids '@:#$r,$a' --new-id-max-allele-len 243 --make-bed --out daphnia.fullids
wc -l daphnia.fullids.bim daphnia.fullids.fam
head -3 daphnia.fullids.bim
/home/go69lef/bin/plink2 --bfile pasteuria --allow-extra-chr --keep daphnia.king354.fam --make-bed --out pasteuria.matched
wc -l pasteuria.matched.bim pasteuria.matched.fam
/home/go69lef/bin/plink2 --bfile pasteuria.matched --allow-extra-chr --geno 0.10 --make-bed --out pasteuria.matched.geno
wc -l pasteuria.matched.geno.bim pasteuria.matched.geno.fam
/home/go69lef/bin/plink2 --bfile pasteuria.matched.geno --allow-extra-chr --maf 0.05 --make-bed --out pasteuria.matched.geno.maf
wc -l pasteuria.matched.geno.maf.bim pasteuria.matched.geno.maf.fam
/home/go69lef/bin/plink2 --bfile pasteuria.matched.geno.maf --allow-extra-chr --set-all-var-ids '@:#$r,$a' --new-id-max-allele-len 243 --make-bed --out pasteuria.fullids
wc -l pasteuria.fullids.bim pasteuria.fullids.fam
head pasteuria.fullids.bim
/home/go69lef/bin/plink2 --bfile pasteuria.fullids --allow-extra-chr --export A --out pasteuria_all_fullids
head -1 pasteuria_all_fullids.raw | cut -f1-12
awk 'BEGIN{OFS="\t"} NR==1 {printf "FID\tIID"; for(i=7;i<=NF;i++) printf "\t%s",$i; printf "\n"; next} {printf "%s\t%s",$1,$2; for(i=7;i<=NF;i++){v=$i; printf "\t%s",(v=="NA" ? "NA" : (v>0 ? 1 : 0))}; printf "\n"}' pasteuria_all_fullids.raw > pasteuria_all_binary.phe
awk 'NR==1 {print "Total columns:",NF; print "Pasteuria phenotypes:",NF-2}' pasteuria_all_binary.phe
wc -l pasteuria_all_binary.phe
head -1 pasteuria_all_binary.phe | tr '\t' '\n' | tail -n +3 | awk '{printf "PAST%06d\t%s\n",NR,$0}' > pasteuria_pheno_map.tsv
head pasteuria_pheno_map.tsv
awk 'BEGIN{OFS="\t"} NR==1 {printf "FID\tIID"; for(i=3;i<=NF;i++) printf "\tPAST%06d",i-2; printf "\n"; next} {print}' pasteuria_all_binary.phe > pasteuria_all_binary_short.phe
head -1 pasteuria_all_binary_short.phe | cut -f1-8
mkdir -p chunks results
n=1
for start in $(seq 3 50 6153); do end=$((start+49)); [ "$end" -gt 6153 ] && end=6153; cut -f"1,2,${start}-${end}" pasteuria_all_binary_short.phe > chunks/chunk_$(printf "%03d" "$n").phe; n=$((n+1)); done
ls chunks/*.phe | wc -l
for i in $(seq -w 1 124); do echo "Running chunk $i"; /home/go69lef/bin/plink2 --bfile daphnia.fullids --allow-extra-chr --pheno chunks/chunk_${i}.phe --1 --covar daphnia.pca.eigenvec --covar-name PC1-PC4 --glm firth-fallback hide-covar skip-invalid-pheno zs --adjust zs --out results/chunk_${i}; rm -f results/chunk_${i}.*.glm.logistic.hybrid.zst; done