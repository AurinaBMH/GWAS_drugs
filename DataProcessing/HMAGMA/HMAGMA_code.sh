WHEREISCODE='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/code/DataProcessing/HMAGMA'
WHEREISGWAS='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/GWAS'
WHEREISOUT='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/GWASgenes'
WHEREIS1000G='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/g1000_eur'
WHEREISANNOT='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/HMAGMA'

# MDD2 dataset: on Chr21 there are 61 SNPs that have only 12 instead of 19 fields. This is because only one cohort
# contributed these SNPs. Remove 61 lines from 21 chromosome,if not removed, will break MAGMA
cat ${WHEREISGWAS}/pgcMDD2.txt | awk '(NF==19)' > ${WHEREISGWAS}/pgcMDD2.txt

for DISORDER in 'ADHD' 'MDD2' 'SCZ' 'DIABETES' 'BIP2' 'HF' 'AD'
do
# with p-values
# 1. use different HMAGMA annotations
for whatANNOT in 'MAGMAdefault' 'Adult_brain' 'Fetal_brain' 'Neuro' 'Astro'
do
${WHEREISCODE}/magma --bfile ${WHEREIS1000G}/g1000_eur --gene-annot ${WHEREISANNOT}/${whatANNOT}.genes.annot --pval ${WHEREISGWAS}/pgc${DISORDER}.txt N=10000 --out ${WHEREISOUT}/pgc${DISORDER}_${whatANNOT}
cp ${WHEREISOUT}/pgc${DISORDER}_${whatANNOT}.genes.out ${WHEREISOUT}/pgc${DISORDER}_${whatANNOT}_genes.txt

# do eMAGMA mapping using eQTL-based annotations send by Zac, these are based on the psychENCODE database
# Options are chosen based on the example on github: https://github.com/eskederks/eMAGMA-tutorial.git

# ###############Gene-based eMAGMA!############
${WHEREISCODE}/magma --bfile ${WHEREIS1000G}/g1000_eur --gene-annot ${WHEREISANNOT}/pec_genes.annot --pval ${WHEREISGWAS}/pgc${DISORDER}.txt N=10000 --gene-settings adap-permp=10000 --out ${WHEREISOUT}/pgc${DISORDER}_eQTLpec
cp ${WHEREISOUT}/pgc${DISORDER}_eQTLpec.genes.out ${WHEREISOUT}/pgc${DISORDER}_eQTLpec_genes.txt
# The --gene-settings flag also controls the settings for permutation-based empirical gene p-values,
# using the fixed-permp or adap-permp modifiers to enable computation of an empirical p-value using
# a fixed number of permutations or an adaptive permutation procedure respectively

done

done

# gene names are in the ENSG00000000457 format, need to get entrez IDs or gene names for them
# use BIOMART annotations to assign IDs to genes in matlab, this is done using get_HMAGMAentrezIDs.m
