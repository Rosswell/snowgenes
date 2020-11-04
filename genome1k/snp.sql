use schema GENOME1K;

-- Staging and permanent tables to handle SNP annotations
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
('NC_000001.10', '1'),
('NC_000002.11', '2'),
('NC_000003.11', '3'),
('NC_000004.11', '4'),
('NC_000005.9', '5'),
('NC_000006.11', '6'),
('NC_000007.13', '7'),
('NC_000008.10', '8'),
('NC_000009.11', '9'),
('NC_000010.10', '10'),
('NC_000011.9', '11'),
('NC_000012.11', '12'),
('NC_000013.10', '13'),
('NC_000014.8', '14'),
('NC_000015.9', '15'),
('NC_000016.9', '16'),
('NC_000017.10', '17'),
('NC_000018.9', '18'),
('NC_000019.9', '19'),
('NC_000020.10', '20'),
('NC_000021.8', '21'),
('NC_000022.10', '22'),
('NC_000023.10', 'X'),
('NC_000024.9', 'Y')
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
