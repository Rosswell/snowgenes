use schema GENOME1K;

-- Table to contain 1000Genomes panel annotation data describing geographic origin and gender for each sample
create or replace table PANEL (
    SAMPLE_ID  varchar,
    POP	        varchar,
    SUPER_POP   varchar,
    GENDER      varchar		
);
