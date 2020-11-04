-- Staging and permanent table to store Gnomad INFO and related fields per-location (excluding Genotype)
use schema GNOMAD;

create or replace table HEADERS (
    VCFNAME varchar,
    LINE    integer,
    TEXT    varchar
);

create or replace table VARIANTS_EXOME_RAW (
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

create or replace table VARIANTS_EXOME (
    CHROM       varchar,
    POS         integer,
    ID          varchar,
    REF         varchar,
    ALT_ALLELE  varchar,
    QUAL        integer,
    FILTER      varchar,
    INFO        variant,
    primary key (CHROM, POS, REF, ALT)
) cluster by (
    CHROM ,
    POS   ,
	REF   
);


create or replace table VARIANTS_GENOME_RAW (
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

create or replace table VARIANTS_GENOME (
    CHROM       varchar,
    POS         integer,
    ID          varchar,
    REF         varchar,
    ALT_ALLELE  varchar,
    QUAL        integer,
    FILTER      varchar,
    INFO        variant,
    primary key (CHROM, POS, REF, ALT)
) cluster by (
    CHROM ,
    POS   ,
	REF   
);

--- Add Schema and tables to the share

GRANT USAGE ON SCHEMA "GENOME"."GNOMAD" TO SHARE "1KGENOMES";
GRANT SELECT ON VIEW "GENOME"."GNOMAD"."HEADERS" TO SHARE "1KGENOMES";
GRANT SELECT ON VIEW "GENOME"."GNOMAD"."VARIANTS_EXOME" TO SHARE "1KGENOMES";
GRANT SELECT ON VIEW "GENOME"."GNOMAD"."VARIANTS_GENOME" TO SHARE "1KGENOMES";
