use schema GENOME1K;

-- Staging and permanent table to store 1000Genomes INFO and related fields per-location (excluding Genotype)
create or replace table VARIANTS_RAW (
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

create or replace table VARIANTS (
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

create or replace SECURE view VARIANTS_FLATTENED as 
select 
    CHROM, 
    POS, 
    ID, 
    REF, 
    ALT, 
    v.index ALT_ALLELE_ID, 
    v.value ALT_ALLELE, 
    info:AC[v.index] AC,
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