# Shell script to preprocess downloaded 1000Genomes data
# Mainly -- break up monolithic large VCF files into smaller subsets, separate headers, and move to S3
# Do this one VCF file at a time to minimize on-disk side from unzipping many VCFs

s3Target='s3://sfc-benchmarks/1000genomes'
#dataDir='/data2/genome'
dataDir='/Users/sozer/desktop/vcf'

# Basic iterator over all chromosome assets
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 \
chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 \
chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrMT;

# For each chromosome, store data in chromosome-specific prefix
do
    echo "processing chromosome ${chr}"
    gzip -d *${chr}.*.gz
    sed -n -f splitvcf.sed *${chr}.*.vcf
    mkdir -p ${chr}
    mv header.txt ${chr}
    mv part*.txt ${chr}
    rm *${chr}.*.vcf
    for file in `ls ${chr}`;
    do
        gzip ${chr}/${file} &
    done
    wait
    aws s3 cp ${dataDir}/${chr} ${s3Target}/${chr} --recursive --include "*" 
    rm -r ${chr}
done

# Handle 1000Genome panel annotation dataset
aws s3 cp ${dataDir}/integrated_call_samples_v3.20130502.ALL.panel  ${s3Target}/

# Handle downloaded NCBI CLINVAR annotation dataset
aws s3 cp ${dataDir}/clinvar.vcf.gz ${s3Target}/

# Handle downloaded NCBI SNP annotation dataset
mv GCF_000001405.25.bgz GCF_000001405.25.gz
bgzip -d GCF_000001405.25.gz
sed -n -f splitsnpvcf.sed GCF_000001405.25
mkdir -p snp
mv header.txt snp
mv part*.txt snp
for file in `ls snp`;
do
    gzip snp/${file} &
done

aws s3 cp ${dataDir}/snp ${s3Target}/snp --recursive --include "*" 