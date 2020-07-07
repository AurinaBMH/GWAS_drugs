WHEREISCODE='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/code/DataProcessing/HMAGMA'
WHEREISGWAS='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/GWAS'
WHEREISCLUMPED='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/GWASclumped'
WHEREISGENE='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/NCBI37.3'
WHEREISOUT='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/GWASgenes'
WHEREIS1000G='/Users/aurinaa/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/g1000_eur'

for DISORDER in 'ADHD' 'MDD2' 'SCZ' 'DIABETES' 'BIP2' 'HF' 'AD'
do
# clump the original GWAS list
echo clumping ${DISORDER} GWAS list
${WHEREISCODE}/plink --noweb --bfile ${WHEREIS1000G}/g1000_eur --clump ${WHEREISGWAS}/pgc${DISORDER}.txt --clump-r2 0.5 --clump-p1 0.00001 --clump-p2 0.01 --out ${WHEREISCLUMPED}/pgc${DISORDER}
# First, do mapping without P-values: take clumped datasets and annotate
# reorder file for MAGMA: SNP, CHR, BP
awk 'BEGIN {FS=" "; OFS="\t"} {print $3, $1, $4}' ${WHEREISCLUMPED}/pgc${DISORDER}.clumped > ${WHEREISCLUMPED}/pgc${DISORDER}_magma.clumped
# run annotation
echo annotating ${DISORDER} GWAS list
${WHEREISCODE}/magma --annotate --snp-loc ${WHEREISCLUMPED}/pgc${DISORDER}_magma.clumped --gene-loc ${WHEREISGENE}/NCBI37.3.gene.loc --out ${WHEREISOUT}/pgc${DISORDER}_genes_clumped
cp ${WHEREISOUT}/pgc${DISORDER}_genes_clumped.genes.annot ${WHEREISOUT}/pgc${DISORDER}_genes_clumped.txt
# pgc${DISORDER}_genes_clumped.txt is the list of genes that map to SNPs based on possition;

# Now, map based on p-values for SNPs, this will give each gene a p-value quantifying the contribution to the phenotype
# 1. make the annotation file based on the whole list of GWAS
# reorder GWASfile for the MAGMA: SNP, CHR, BP
if DISORDER='AD'
then
  # the order of columns in the AD file is different from others
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $6, $2, $3}' ${WHEREISGWAS}/pgc${DISORDER}.txt > ${WHEREISGWAS}/pgc${DISORDER}_magma.txt
else
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $2, $1, $3}' ${WHEREISGWAS}/pgc${DISORDER}.txt > ${WHEREISGWAS}/pgc${DISORDER}_magma.txt
fi
#echo making gene file for ${DISORDER} GWAS list
${WHEREISCODE}/magma --annotate --snp-loc ${WHEREISGWAS}/pgc${DISORDER}_magma.txt --gene-loc ${WHEREISGENE}/NCBI37.3.gene.loc --out ${WHEREISOUT}/pgc${DISORDER}_genes_all
# with p-values
${WHEREISCODE}/magma --bfile ${WHEREIS1000G}/g1000_eur --gene-annot ${WHEREISANNOT}/pgc${DISORDER}_genes_all.genes.annot --pval ${WHEREISGWAS}/pgc${DISORDER}.txt N=10000 --out ${WHEREISOUT}/pgc${DISORDER}
cp ${WHEREISOUT}/pgc${DISORDER}.genes.out ${WHEREISOUT}/pgc${DISORDER}_genes_all.txt

done
