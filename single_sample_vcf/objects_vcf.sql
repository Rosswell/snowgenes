-- Basic components of a schema to ingest VCF data containing a SINGLE SAMPLE PER FILE.
-- Each row contains the FILENAME for lineage and to potentially link back to sample metadata.
-- This schema can be amended to extract and maintain a SAMPLE_ID if it is parseable from
-- the filename, or from information in headers if needed. 

use schema VCF;

-- Stage for VCF files
create or replace stage VCF_STAGE
    URL = 's3://**** BUCKET PATH ****'
    CREDENTIALS = (  AWS_KEY_ID =  '*****' AWS_SECRET_KEY = '*****' )
    FILE_FORMAT = (TYPE = CSV FIELD_DELIMITER = '\t' COMPRESSION = AUTO ERROR_ON_COLUMN_COUNT_MISMATCH = FALSE);

-- format for vcf headers
 create or replace file format vcf_stage_fmt  TYPE = CSV COMPRESSION = AUTO FIELD_DELIMITER = '`';
 
 -- Raw staging table to contain all rows + headers from a gvcf
 -- Should be truncated after each batch load
create or replace table VCF_RAW_STG (
    filename varchar,
    line     integer,
    text     varchar
);

create or replace table HEADERS (
    filename  varchar,
    line      integer,
    text      varchar
);

create or replace table GENOTYPES (
    CHROM       varchar,
    POS         integer,
    ID          varchar,
    REF         varchar,
    ALT         array  , 
    QUAL        integer,
    FILTER      varchar,
    INFO        variant, 
    FORMAT      varchar,
    GT          varchar(4),
    FILENAME    varchar,
    primary key (CHROM, POS)
) cluster by (
    CHROM ,
    POS   
 );
