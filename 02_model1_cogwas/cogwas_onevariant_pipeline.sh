#!/usr/bin/env bash
set -euo pipefail
cd /Users/simalardil/Desktop/coGWAS_paperstyle
pwd
/Users/simalardil/plink2 --version
ls -lh *.vcf* 
file daphnia_input.vcf.gz
/Users/simalardil/plink2 --vcf daphnia_input.vcf.gz --double-id --allow-extra-chr --max-alleles 2 --make-bed --out daphnia
 wc -l daphnia.bim daphnia.fam
/Users/simalardil/plink2 --vcf pasteuria_input.vcf.gz --double-id --allow-extra-chr --max-alleles 2 --make-bed --out pasteuria
wc -l pasteuria.bim pasteuria.fam
/Users/simalardil/plink2 --bfile daphnia --allow-extra-chr --geno 0.10 --make-bed --out daphnia.geno
wc -l daphnia.geno.bim daphnia.geno.fam
/Users/simalardil/plink2 --bfile daphnia.geno --allow-extra-chr --mind 0.10 --make-bed --out daphnia.geno.mind
wc -l daphnia.geno.mind.bim daphnia.geno.mind.fam
/Users/simalardil/plink2 --bfile daphnia.geno.mind --allow-extra-chr --maf 0.05 --make-bed --out daphnia.geno.mind.maf
wc -l daphnia.geno.mind.maf.bim daphnia.geno.mind.maf.fam
/Users/simalardil/plink2 --bfile daphnia.geno.mind.maf --allow-extra-chr --king-cutoff 0.354 --make-bed --out daphnia.king354
wc -l daphnia.king354.bim daphnia.king354.fam daphnia.king354.king.cutoff.out.id
/Users/simalardil/plink2 --bfile daphnia.king354 --allow-extra-chr --set-all-var-ids '@:#$r,$a' --new-id-max-allele-len 243 --indep-pairwise 200kb 1 0.2 --out daphnia.prune
wc -l daphnia.prune.prune.in daphnia.prune.prune.out
/Users/simalardil/plink2 --bfile daphnia.king354 --allow-extra-chr --set-all-var-ids '@:#$r,$a' --new-id-max-allele-len 243 --extract daphnia.prune.prune.in --pca 10 --out daphnia.pca
/Users/simalardil/plink2 --bfile pasteuria --allow-extra-chr --keep daphnia.king354.fam --make-bed --out pasteuria.matched
wc -l pasteuria.matched.bim pasteuria.matched.fam
/Users/simalardil/plink2 --bfile pasteuria.matched --allow-extra-chr --geno 0.10 --make-bed --out pasteuria.matched.geno
wc -l pasteuria.matched.geno.bim pasteuria.matched.geno.fam
/Users/simalardil/plink2 --bfile pasteuria.matched.geno --allow-extra-chr --maf 0.05 --make-bed --out pasteuria.matched.geno.maf
wc -l pasteuria.matched.geno.maf.bim pasteuria.matched.geno.maf.fam
/Users/simalardil/plink2 --bfile pasteuria --allow-extra-chr --chr 1 --from-bp 4530 --to-bp 4530 --export A --out pasteuria_chr1_4530_before
awk 'NR>1 {count[$7]++} END {for (value in count) print value, count[value]}' pasteuria_chr1_4530_before.raw | sort -n
/Users/simalardil/plink2 --bfile pasteuria.matched.geno.maf --allow-extra-chr --chr 1 --from-bp 4530 --to-bp 4530 --export A --out pasteuria_chr1_4530_after
awk 'NR>1 {count[$7]++} END {for (value in count) print value, count[value]}' pasteuria_chr1_4530_after.raw | sort -n
awk 'BEGIN{OFS="\t"; print "FID","IID","PASTEURIA_G_PRESENT"} NR>1 {v=$7; print $1,$2,(v=="NA" ? "NA" : (v>0 ? 1 : 0))}' pasteuria_chr1_4530_after.raw > pheno.txt
/Users/simalardil/plink2 --bfile daphnia.king354 --allow-extra-chr --pheno pheno.txt --pheno-name PASTEURIA_G_PRESENT --1 --covar daphnia.pca.eigenvec --covar-name PC1-PC4 --glm firth-fallback hide-covar --adjust --out cogwas_onevariant
