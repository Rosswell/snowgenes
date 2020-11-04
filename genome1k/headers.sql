use schema GENOME1K;

-- Table to store all header information from loaded VCFs

create or replace table HEADERS (
    VCFNAME varchar,
    LINE    integer,
    TEXT    varchar
);
