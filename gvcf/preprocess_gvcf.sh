#!/bin/bash
# Convert BCF files into compressed gVCF header and detail files for loading into Snowflake
# This alternative load pipeline is designed to be faster, but requires preprocessing
#
# This demonstration code assumes samples are named   1.bcf - 100.bcf
# 
# Preprocessing locates the SAMPLE_ID in the last line of each VCF Header, and prepends the
#   ID to the filename so that it can be used at load time to annotate each row
#
# Preprocessing should be followed with executing load_gvcf.sql
#
s3Target='s3://*****'
#dataDir='/data2/genome'
dataDir='/Users/*****'

for sample in {1..100};
do
    bcftools convert -O v ${sample}.bcf | grep -v '#' | gzip > ${sample}.vcfrows.gz
    bcftools convert -O v ${sample}.bcf | head -5000 | grep  '#' | gzip  > ${sample}.vcfheaders.gz
    newname=`gzip -d -c ${sample}.vcfheaders.gz | grep '#CHROM.*' | cut -f 10`
    mv ${sample}.vcfrows.gz ${newname}.${sample}.vcfrows.gz
    mv ${sample}.vcfheaders.gz ${newname}.${sample}.vcfheaders.gz
done
mkdir -p preprocessed
mv *.gz preprocessed

aws s3 cp ${datadir}/preprocessed ${s3Target}/preprocessed --recursive --exclude "*" --include "*.gz"

# Handle downloaded NCBI SNP annotation dataset
mv GCF_000001405.38.bgz GCF_000001405.38.gz
bgzip -d GCF_000001405.38.gz
sed -n -f splitsnpvcf.sed GCF_000001405.38
mkdir -p snp
mv header.txt snp
mv part*.txt snp
for file in `ls snp`;
do
    gzip snp/${file} &
done

aws s3 cp ${dataDir}/snp ${s3Target}/snp --recursive --include "*" 

# Handle Clinvar
aws s3 cp ${dataDir}/clinvar.vcf.gz ${s3Target}/