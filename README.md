# SnowGenes

v2 August 2019 – Stuart Ozer – Snowflake Computing
#

## Overview

SnowGenes is a demonstration project and data sharing example designed to showcase Snowflake's capabilities in working with genomic analytical workflows and data sources.  Our initial examples consider 2 types of sources:

- 1000Genome Project Data and related annotations
- gVCF data from multiple individuals

These are described in more detail below.

#
## 1000Genome Project Data
 
We have created a sample Genomics schema to demonstrate how Snowflake can be used to manage information generated by modern genomic sequencing and annotation technologies.  The Multi-TB scale dataset is designed to support high-performance querying for patterns that are common in analytical and data exploration workflows from the genomics community.
 
The data includes:
 
1.	**1000Genomes Phase 3 Full Data**.  This is the complete contents of Phase 3 sequencing from the 1000 Genomes project, representing over 2500 distinct individuals, including all Chromosomes and Mitochondrial DNA from the subjects.   

    This corresponds to the raw data found in the VCF files at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/.  Most raw data is dated 5/27/15.  It utilizes position numbers and reference allele values based on  reference genome GRCh37, also known as hg19.
 
    The 1000Genome data includes some useful annotations including the frequencies of individual variations found in different populations, and the geographic origin of each of the sample genomes.   
 
2.	**NCBI Clinvar database**.   A wide set of annotations including disease and gene associations relevant to each variant.  The reference genome and positions are based on GRCh37, and were downloaded from  https://www.ncbi.nlm.nih.gov/clinvar/
 
3.	**NCBI SNP database**.   Documentation and references annotating most SNPs found in the 1000Genomes.  Again, the reference genome and positions are based on GRCh37, and the dataset is downloadable as GCF_000001405.25.bgz at ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF 
 
Information from the 3 datasets can be used together since they have common reference genome basis.
 
## Schema
 
### 1000Genomes dataset
 
Rather than storing raw VCF data from 1000Genomes, we have performed a series of light transformations on the VCF that makes querying more efficient and natural in a SQL-based environment.  
 
#### Genotypes
 
The most important innovation we have made for the large genome dataset is that we provide a unique row in the database for every individual, at every Variant position in each chromosome.  Basically, we "UNPIVOT" the raw 1000Genomes VCF data.  But we also order (cluster) this data by location (chromosome & position) and by genotype.  As a result the final dataset is relatively compact in Snowflake's columnar store, while allowing for very efficient location-specific or range-based queries.  This data is contained in the large Genotypes table, taking up about 1.3 TB compressed space.  This translates to approx 5 bytes per row inside Snowflake, or ~500MB per genome.

We also materialize the values for Allele1 and Allele2 at each position, decoding the GT field to minimize computational overhead during large population queries.
 
This organization of Genotypes is especially efficient for seeking individuals or evaluating statistics on variants in a specific location or range of locations within the genome.  But it is not efficient for testing a single individual's genome for indicators or features genome-wide.   We take advantage of Snowflake's Materialized View feature to maintain a second 'copy' of Genotypes clustered by Individual.  This view is named Genotypes_By_Sample.
```
create or replace table GENOTYPES (
   CHROM varchar,
   POS integer,
   ID array ,
   REF varchar,
   ALT array ,
   SAMPLE_ID varchar,
   GT varchar,
   ALLELE1 varchar,
   ALLELE2 varchar
primary key (CHROM, POS, REF, ALT, SAMPLE_ID)
) cluster by (
   CHROM ,
   POS ,
   REF ,
   GT
);
 
 
create or replace materialized view GENOTYPES_BY_SAMPLE (
   CHROM, 
   POS, 
   ID, 
   REF, 
   ALT, 
   SAMPLE_ID, 
   GT,
   ALLELE1,
   ALLELE2
)
cluster by (
   SAMPLE_ID, 
   CHROM,
   POS,
   REF
)
as select CHROM, POS, ID, REF, ALT, SAMPLE_ID, GT, ALLELE1, ALLELE2
from genotypes
;
``` 
 
#### Variants 
 
We also separate the informational annotations about locations (contained in the INFO field) -- such as allele frequencies, population-specific occurrence rates, etc. -- into a separate table with one entry per variant position.  This is the Variants table, and it can be easily joined to Genomes based on the common key (CHROM, POS, REF, ALT).  Variants also contains the light annotation-per-position found in the REF and FILTER fields of the original 1000Genomes VCF.
 
Much of the information in Variants is stored in arrays within the INFO field -- with one entry for each comma-separated variant present in the ALT field.  Examples are population frequencies per ALT-contained value.  For much analysis it makes sense to "Flatten" this arrayed data, obtaining one-row-per-variant-allele.  To simplify this common pattern for many queries, we create a (non-materialized) view Variants_Flattened that exposes a single allele in each row, as well as the allele-specific values from the INFO field.
```
create or replace table VARIANTS (
   CHROM varchar,
   POS integer,
   ID array ,
   REF varchar,
   ALT array ,
   QUAL integer,
   FILTER varchar,
   INFO variant,
primary key (CHROM, POS, REF, ALT)
) cluster by (
   CHROM ,
   POS ,
   REF 
);
 
 
create or replace view VARIANTS_FLATTENED as 
select 
   CHROM, 
   POS, 
   REF, 
   ALT, 
   v.index ALLELE_ID, 
   v.value ALLELE, 
   info:AC[v.index] AC
   info:AF[v.index] AF, 
   info:EAS_AF[v.index] EAS_AF, 
   info:AMR_AF[v.index] AMR_AF, 
   info:AFR_AF[v.index] AFR_AF, 
   info:EUR_AF[v.index] EUR_AF, 
   info:SAS_AF[v.index] SAS_AF,
   info INFO
from variants,
lateral flatten (INPUT => alt, mode => 'ARRAY') v
;

```
 
### Annotation data 
 
#### DbSNP
 
Data about SNPs, downloaded from NCBI, are stored in the table named SNP.   This data can be joined to Variants and/or Genotypes on columns CHROM, POS, REF -- which serve as an approximate key.  In SNP, there may be an exact match, subset, superset, or no match of ALT values at a position compared to the ALT value in 1000Genomes Variants and Genomes table -- so to find relevant matching rows it is best to add 

`WHERE arrays_overlap (SNP.ALT, VARIANTS.ALT)` 

or similar comparison to the join condition.
``` 
create or replace table SNP (
   CHROM varchar,
   ORIG_CHROM varchar,
   POS integer,
   ID array ,
   REF varchar,
   ALT array ,
   QUAL integer,
   FILTER varchar,
   INFO variant,
primary key (CHROM, POS, REF, ALT)
) cluster by (
   CHROM ,
   POS ,
   REF 
);
```
 
#### Clinvar
 
The assembled information about clinical indicators, also downloaded from NCBI, is stored in table Clinvar.  Like SNP, it can be joined to Variants and Genotypes on the common key CHROM, POS, REF.  Unlike both SNP and Variants (and Genotypes), the Clinvar table specifies only a single ALT value per row.  If different ALT alleles are annotated at the same position, there are other rows to cover them.   So to find relevant matching rows for Variants or Genotypes, add a condition such as 
 
`WHERE arrays_overlap (CLINVAR.ALT, VARIANTS.ALT)`
``` 
create or replace table CLINVAR (
   CHROM varchar,
   POS integer,
   ID array ,
   REF varchar,
   ALT array ,
   QUAL integer,
   FILTER varchar,
   INFO variant,
primary key (CHROM, POS, REF, ALT)
) cluster by (
   CHROM ,
   POS ,
   REF 
);
```
 
### VCF Headers
 
#### Headers 
 
VCF files all begin with header rows that define, among other things, the possible fields contained inside the INFO column of the associated data.  These headers differ for each data source loaded from VCF and serve as a sort of dynamic documentation for the data rows in other tables.
 
Our schema maintains all header information from loaded VCF tables (whether from genomes or annotation sources) in a table called Headers.   Info from Headers can be interpreted using documentation such as https://samtools.github.io/hts-specs/VCFv4.2.pdf.
 
Headers from all VCF sources are contained in the same table.  A specific data source is identified by the column vcfname.  In our schema, possible values for vcfname include:
```
clinvar
snp
chr1
chr2
.
.
.
chr21
chrX
chrY
chrMT
```
The various chr… headers describe info field contents for the Variants table which was populated by the data from 1000Genomes individual chromosomes' VCFs.
 
For a specific value of vcfname, header line numbers are indicated by the column line.  So 
 
`Select * from headers where vcfname='snp' order by line`
 
Will return the full set of header lines from the DbSNP VCF file header. 
```
create or replace table HEADERS (
   VCFNAME varchar,
   LINE integer,
   TEXT varchar
);
```
 
## Helper Functions
 
We have included a number of User Defined Functions (UDFs) to make it simpler to parse and extract data from tables derived from VCFs.
 
These UDFs include:
 
`Allele1(ref varchar, alt variant, gt varchar)`
 
Returns the character value of the first allele in a Genotype GT field.  For example, for a Genotypes row with ref='C', alt= ['T'] and gt = '1|0', `Allele1(ref, alt, gt) = 'T'`
 
`Allele2(ref varchar, alt variant, gt varchar)`
 
Similar to Allele1, but delivers the second allele.  In a haploid sample, it will return NULL.   For a Genotypes row with ref='C', alt= ['T'] and gt = '1|0', `Allele2(ref, alt, gt) = 'C'`
 
`allele_position(gt varchar, which int)`
 
To support Array Lookups into INFO fields containing separate entries per ALT value, this function takes a genotype and the strand number (0 or 1) and returns the integer position in the ALT array or any other allele-specific array corresponding to that allele.  If the allele is REF valued (e.g. noted as 0), the function returns NULL.   For a Genotypes row with ref='C', alt= ['T', 'G'] and gt = '2|0', `Allele_position(gt, 0) = 1` and
`Allele_position(gt, 1) = null`
 
Actual values for the 2 alleles, or relevent values for any other array that is allele-specific, can be retrieved using these position indexes.  For example
```
select *, 
coalesce(alt[allele_position(gt,0)], ref), coalesce(alt[allele_position(gt,1)], ref) 
from genotypes …
```
returns the character values of the 2 strands at any position
 
`is_homozygous (gt varchar)`
 
Boolean that returns TRUE for a genotype value where both strands have same value at a position (or is single-stranded), or FALSE otherwise.
```
is_homozygous('1|1') = TRUE
is_homozygous('0|1') = FALSE
is_homozygous('0')   = TRUE
```

#
## gVCF Data

**coming soon ...** 
