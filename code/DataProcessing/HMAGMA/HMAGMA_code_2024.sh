WHEREISCODE='/Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/code/DataProcessing/HMAGMA'
WHEREISGWAS='/Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/data/GWASlists/GWAS_2022'
WHEREISOUT='/Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/data/GWASlists/GWASgenes_2022'
WHEREIS1000G='/Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/data/GWASlists/g1000_eur'
WHEREISANNOT='/Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/data/GWASlists/HMAGMA'

# SCZ3 GWAS contains a line with missing values, it needed to be removed for MAGMA to run
# line: 21 rs148878475 9648204 C T 0.9878 0.9845 0.1375 7.5732 1.0014 0.0432 0
# MDD2 dataset: on Chr21 there are 61 SNPs that have only 12 instead of 19 fields. This is because only one cohort
# contributed these SNPs. Remove 61 lines from 21 chromosome,if not removed, will break MAGMA
# cat ${WHEREISGWAS}/pgcMDD2.txt | awk '(NF==19)' > ${WHEREISGWAS}/pgcMDD2.txt
# DIABETES GWAS contains a line with a very low p-value, that is read as non-number;
# it needed to be removed for MAGMA to run
# line: 10 114758349 rs7903146 T C 0.291584512494863 0.3059 0.0077 1.33e-347 597475

# DIABETES2 GWAS contains a few lines with a very low p-value, that is read as non-number; Should be 1e-300, but filtering automatically
# also doesn't read such values as numbers. Therefore, this is done manually. These lines needed to be removed for MAGMA to run
# remove one by one;

# sed -i '' '6660720d' pgcDIABETES2.txt
# sed -i '' '6660734d' pgcDIABETES2.txt
# sed -i '' '6660735d' pgcDIABETES2.txt
# sed -i '' '6660736d' pgcDIABETES2.txt
# sed -i '' '6660744d' pgcDIABETES2.txt
# sed -i '' '6660741d' pgcDIABETES2.txt
# sed -i '' '6660761d' pgcDIABETES2.txt
# sed -i '' '6660763d' pgcDIABETES2.txt
# sed -i '' '6660774d' pgcDIABETES2.txt
# sed -i '' '6660824d' pgcDIABETES2.txt
# sed -i '' '6660832d' pgcDIABETES2.txt
# sed -i '' '6660864d' pgcDIABETES2.txt
# sed -i '' '6660867d' pgcDIABETES2.txt
# sed -i '' '6660896d' pgcDIABETES2.txt
# sed -i '' '6660920d' pgcDIABETES2.txt


# latest dataset main text: 'MDD4' 'ADHD3' 'SCZ3' 'BIP3' 'DIABETES2' 'HF' 'RA' 'IBD'
# older datasets 2021 option: 'MDD3' 'ADHD' 'SCZ' 'BIP2' 'DIABETES'
for DISORDER in 'MDD4' 'ADHD3' 'SCZ3' 'BIP3' 'DIABETES2' 'HF' 'RA' 'IBD' 'MDD3' 'ADHD' 'SCZ' 'BIP2' 'DIABETES'
do
# with p-values
# 1. use different HMAGMA annotations
# 'MAGMAdefault' 'Adult_brain' 'Fetal_brain' 'Neuro' 'Astro'
for whatANNOT in 'MAGMAdefault' 'Adult_brain' 'Fetal_brain' 'Neuro' 'Astro'
do
${WHEREISCODE}/magma --bfile ${WHEREIS1000G}/g1000_eur --gene-annot ${WHEREISANNOT}/${whatANNOT}.genes.annot --pval ${WHEREISGWAS}/pgc${DISORDER}.txt N=10000 --out ${WHEREISOUT}/pgc${DISORDER}_${whatANNOT}
cp ${WHEREISOUT}/pgc${DISORDER}_${whatANNOT}.genes.out ${WHEREISOUT}/pgc${DISORDER}_${whatANNOT}_genes.txt
done

# do eMAGMA mapping using eQTL-based annotations send by Zac, these are based on the psychENCODE database
# Options are chosen based on the example on github: https://github.com/eskederks/eMAGMA-tutorial.git
# eMAGMA using brain annotations

${WHEREISCODE}/magma --bfile ${WHEREIS1000G}/g1000_eur --gene-annot ${WHEREISANNOT}/pec_genes.annot --pval ${WHEREISGWAS}/pgc${DISORDER}.txt N=10000 --gene-settings adap-permp=10000 --out ${WHEREISOUT}/pgc${DISORDER}_eQTLbrain
cp ${WHEREISOUT}/pgc${DISORDER}_eQTLbrain.genes.out ${WHEREISOUT}/pgc${DISORDER}_eQTLbrain_genes.txt

# run separately for eMAGMA using blood, heart, liver eQTLs (from GTEx)
# 'Small_Intestine_Terminal_Ileum' 'Pancreas' 'Whole_Blood' 'Liver' 'Heart_Left_Ventricle' 'Colon_Transverse' 'Colon_Sigmoid' 'Adipose_Subcutaneous' 'Adipose_Visceral_Omentum'
for whatANNOT in 'Small_Intestine_Terminal_Ileum' 'Pancreas' 'Whole_Blood' 'Liver' 'Heart_Left_Ventricle' 'Colon_Transverse' 'Colon_Sigmoid' 'Adipose_Subcutaneous' 'Adipose_Visceral_Omentum'
do
# ###############Gene-based eMAGMA!############
${WHEREISCODE}/magma --bfile ${WHEREIS1000G}/g1000_eur --gene-annot ${WHEREISANNOT}/emagma_gtexv8_annot/${whatANNOT}.genes.annot --pval ${WHEREISGWAS}/pgc${DISORDER}.txt N=10000 --gene-settings adap-permp=10000 --out ${WHEREISOUT}/pgc${DISORDER}_eQTL${whatANNOT}
cp ${WHEREISOUT}/pgc${DISORDER}_eQTL${whatANNOT}.genes.out ${WHEREISOUT}/pgc${DISORDER}_eQTL${whatANNOT}_genes.txt

done


done

# gene names are in the ENSG00000000457 format, need to get entrez IDs or gene names for them
# use BIOMART annotations to assign IDs to genes in matlab, this is done using get_HMAGMAentrezIDs.m
