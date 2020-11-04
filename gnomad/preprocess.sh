# Shell script to both download and preprocess GNOMAD data
# Mainly -- break up monolithic large VCF files into smaller subsets, separate headers, and move to S3
# Do this one VCF file at a time to minimize on-disk side from unzipping many VCFs
# For storage space efficiency, Download each chromosome and process, delete before loading next chromosome

# Perform separate iterations for the EXOME and GENOME data sets in GNOMAD
for dataset in exome genome;
do

    s3Target="s3://sfc-benchmarks/gnomad/${dataset}"
    dataDir='/data2/genome'

    # For each exome chromosome, store data in chromosome-specific prefix
    # Basic iterator over all chromosome assets
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 \
    chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 \
    chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;
    do
        echo "processing chromosome ${chr}"
        wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/${dataset}s/gnomad.${dataset}s.r2.1.sites.${chr}.vcf.bgz
        mv gnomad.${dataset}s.r2.1.sites.${chr}.vcf.bgz gnomad.${dataset}s.r2.1.sites.${chr}.vcf.gz
        gzip -d *${chr}.*.gz
        sed -n -f split${dataset}vcf.sed *${dataset}s*${chr}.vcf
        mkdir -p ${chr}
        mv header.txt ${chr}
        mv part*.txt ${chr}
        rm *${dataset}s*${chr}.vcf
        for file in `ls ${chr}`;
        do
            gzip ${chr}/${file} &
        done
        wait
        aws s3 cp ${dataDir}/${chr} ${s3Target}/${chr} --recursive --include "*" 
        rm -r ${chr}
    done
done