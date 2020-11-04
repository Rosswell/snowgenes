-- Tables to hold genotype data with one row per 
-- individual per position.  Clustered by position to support 
-- efficient lookups and statistics-gathering per-position.
--
-- Included is a materialized view for the Sample-specific alternative clustering order
-- supporting efficient lookups by sample_id.

use schema GENOME1K;


create or replace table GENOTYPES_STG (
    CHROM       char(2),
    POS         integer,
    ID          array  ,
    REF         varchar(256),
    ALT         array  ,
    SAMPLE_ID   varchar(8),
    GT          varchar(8),
    primary key (CHROM, POS, REF, ALT, SAMPLE_ID)
) cluster by (
    CHROM ,
    POS   ,
	REF   ,
    GT
 );



-- ************** I M P O R T A N T ******************************************
-- DEFER EXECUTION OF THESE STATEMENTS UNTIL GENOTYPES_STG TABLE IS POPULATED! 
-- E.g. execute this snippet after completing loading.sql
-- ***************************************************************************
create or replace table GENOTYPES (
    CHROM       char(2),
    POS         integer,
    ID          array  ,
    REF         varchar(256),
    ALT         array  ,
    SAMPLE_ID   varchar(8),
    GT          varchar(8),
    ALLELE1     varchar,
    ALLELE2     varchar,
    primary key (CHROM, POS, REF, ALT, SAMPLE_ID)
) cluster by (
    CHROM ,
    POS   ,
	REF   ,
    GT
 ) as
    SELECT CHROM, POS, ID, REF, ALT, SAMPLE_ID, GT,
       ALLELE1(ref, alt, gt) as ALLELE1,
       ALLELE2(ref, alt, gt) as ALLELE2
    FROM GENOTYPES_STG
;

-- It is much faster and cheaper to build the MV using populated source data
-- rather than using lazy creation as base table is populated
create or replace secure materialized view GENOTYPES_BY_SAMPLE(
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
from GENOTYPES
;
-- *****************************************************************************