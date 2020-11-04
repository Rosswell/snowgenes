-- ALTERNATIVE GVCF schema designed around higher-performance loading
-- Requires some preprocessing of input file content prior to load
-- Leaves more semi-structured data un-parsed until query time
--
use schema GVCF;

create or replace table HEADERS (
    FILENAME  varchar,
    SAMPLE_ID varchar,
    LINE      integer,
    TEXT      varchar
);


create or replace table GENOTYPES_RAW (
    CHROM       char(2),
    POS         integer,
    END         integer,
    ID          varchar(64),
    REF         varchar(256),
    ALT_STRING  varchar, 
    QUAL        integer,
    FILTER      varchar(32),
    INFO        varchar,
    SAMPLE_ID   varchar(11),
    FORMAT      varchar(128),
    GT          varchar(4),
    VALS        varchar,
    FILENAME    varchar(22)
) 
;

create or replace view GENOTYPES_BY_SAMPLE as 
select 
    CHROM as CHROM, 
    POS as POS, 
    END as END,
    ID as ID,
    REF as REF,
    ALT_STRING as ALT_STRING,
    split(ALT_STRING, ',') as ALT,
    QUAL as QUAL,
    FILTER as FILTER,
    SAMPLE_ID as SAMPLE_ID,
    FORMAT as FORMAT,
    INFO as INFO,
    GT as gt,
    VALS as VALS,
    strtok(VALS, ':', 1 + array_position('DP'::variant, split(FORMAT, ':')))::integer as DP,
    strtok(VALS, ':', 1 + array_position('GQ'::variant, split(FORMAT, ':')))::integer as GQ,
    strtok(VALS, ':', 1 + array_position('MIN_DP'::variant, split(FORMAT, ':')))::integer as MIN_DP,
    strtok(VALS, ':', 1 + array_position('PS'::variant, split(FORMAT, ':')))::integer as PS,
    strtok(VALS, ':', 1 + array_position('VAR_TYPE'::variant, split(FORMAT, ':')))::varchar as VAR_TYPE,
    strtok(VALS, ':', 1 + array_position('VAR_CONTEXT'::variant, split(FORMAT, ':')))::varchar as VAR_CONTEXT,
    strtok(VALS, ':', 1 + array_position('STR_MAX_LEN'::variant, split(FORMAT, ':')))::integer as STR_MAX_LEN,
    strtok(VALS, ':', 1 + array_position('STR_PERIOD'::variant, split(FORMAT, ':')))::integer as STR_PERIOD,
    strtok(VALS, ':', 1 + array_position('STR_TIMES'::variant, split(FORMAT, ':')))::float as STR_TIMES,
    FILENAME as FILENAME
from genotypes_raw
;

create or replace materialized view GENOTYPES_BY_POSITION (
    CHROM       ,
    POS         ,
    END         ,
    ID          ,
    REF         ,
    ALT_STRING  ,
    ALT         ,
    QUAL        ,
    FILTER      ,
    SAMPLE_ID   ,
    FORMAT      ,
    INFO        ,
    GT          ,
    DP          ,
    GQ          ,
    MIN_DP      ,
    PS          ,
    VAR_TYPE    ,
    VAR_CONTEXT ,
    STR_MAX_LEN ,
    STR_PERIOD  ,
    STR_TIMES   ,
    FILENAME    
)
cluster by (
    CHROM,
    POS
 ) as 
select 
    CHROM as CHROM, 
    POS as POS, 
    END as END,
    ID as ID,
    REF as REF,
    ALT_STRING as ALT_STRING,
    split(ALT_STRING, ',') as ALT,
    QUAL as QUAL,
    FILTER as FILTER,
    SAMPLE_ID as SAMPLE_ID,
    FORMAT as FORMAT,
    INFO as INFO,
    GT as gt,
    strtok(VALS, ':', 1 + array_position('DP'::variant, split(FORMAT, ':')))::integer as DP,
    strtok(VALS, ':', 1 + array_position('GQ'::variant, split(FORMAT, ':')))::integer as GQ,
    strtok(VALS, ':', 1 + array_position('MIN_DP'::variant, split(FORMAT, ':')))::integer as MIN_DP,
    strtok(VALS, ':', 1 + array_position('PS'::variant, split(FORMAT, ':')))::integer as PS,
    strtok(VALS, ':', 1 + array_position('VAR_TYPE'::variant, split(FORMAT, ':')))::varchar as VAR_TYPE,
    strtok(VALS, ':', 1 + array_position('VAR_CONTEXT'::variant, split(FORMAT, ':')))::varchar as VAR_CONTEXT,
    strtok(VALS, ':', 1 + array_position('STR_MAX_LEN'::variant, split(FORMAT, ':')))::integer as STR_MAX_LEN,
    strtok(VALS, ':', 1 + array_position('STR_PERIOD'::variant, split(FORMAT, ':')))::integer as STR_PERIOD,
    strtok(VALS, ':', 1 + array_position('STR_TIMES'::variant, split(FORMAT, ':')))::float as STR_TIMES,
    FILENAME as FILENAME
from genotypes_raw
;

 -- Staging and permanent tables to support Clinvar annotation in HG38 format
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

-- Staging and permanent tables to handle SNP annotations in HG38 format
create or replace table SNP_RAW (
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

-- Translation trable from Genbank chromosome IDs used in SNP DB, to standard human Chrom #s.
create or replace table GENBANK_TO_CHROM (
    GENBANK_ID varchar, 
    CHROM varchar);

insert into GENBANK_TO_CHROM (GENBANK_ID, CHROM) values 
('NC_000001.11', '1'),
('NC_000002.12', '2'),
('NC_000003.12', '3'),
('NC_000004.12', '4'),
('NC_000005.10', '5'),
('NC_000006.12', '6'),
('NC_000007.14', '7'),
('NC_000008.11', '8'),
('NC_000009.12', '9'),
('NC_000010.11', '10'),
('NC_000011.10', '11'),
('NC_000012.12', '12'),
('NC_000013.11', '13'),
('NC_000014.9', '14'),
('NC_000015.10', '15'),
('NC_000016.10', '16'),
('NC_000017.11', '17'),
('NC_000018.10', '18'),
('NC_000019.10', '19'),
('NC_000020.11', '20'),
('NC_000021.9', '21'),
('NC_000022.11', '22'),
('NC_000023.11', 'X'),
('NC_000024.10', 'Y')
;

create or replace table SNP (
    CHROM       varchar,
    ORIG_CHROM  varchar,
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

-- Example loading of BED file and conversion into BED_POSITIONS table (one row per variant position of interest)
create or replace table exac_bed (
    CHROM       varchar(2),
    POS         integer,
    END_PLUS1   integer
) cluster by (
    CHROM,
    POS
);

