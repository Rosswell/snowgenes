-- Staging and Permanent tables to hold Clinvar VCF files 
-- Raw data downloadable from https://www.ncbi.nlm.nih.gov/clinvar/
use schema GENOME1K;

create or replace table CLINVAR_RAW (
    CHROM       varchar,
    POS         integer,
    ID          varchar,
    REF         varchar,
    ALT         varchar,
    QUAL        integer,
    FILTER      varchar,
    INFO        varchar,
    primary key (CHROM, POS, REF, ALT)
) cluster by (
    CHROM ,
    POS   ,
	REF   ,
	ALT   
 );

create or replace table CLINVAR (
    CHROM       varchar,
    POS         integer,
    ID          array  ,
    REF         varchar,
    ALT         array  ,
    QUAL        integer,
    FILTER      varchar,
    INFO        variant,
    primary key (CHROM, POS, REF, ALT)
) cluster by (
    CHROM ,
    POS   ,
	REF   
 );
