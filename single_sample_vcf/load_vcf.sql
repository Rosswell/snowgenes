-- Load script to populate GENOTYPES and HEADERS from a batch of VCF files
-- in a stage.  Loading uses a staging file to ingest data prior to parsing 
-- INFO fields for the final load.
--   
-- In this version, only GT is assumed to be present in the FORMAT field of VCF data.   
-- Schema and script can be amended to add fully structured or additional semi-structured
-- columns to support more complex FORMAT entries.

use schema VCF;

-- Remove all rows from the staging table
-- We use DELETE instead of TRUNCATE to avoid loading staging files 
-- multiple times as they accumulate.  (TRUNCATE resets tracking for previously 
-- loaded files forcing administrator to remove loaded files from S3 prior to loading a new batch 
delete from vcf_raw_stg;

-- Load a staging version of the vcf 
copy into vcf_raw_stg from  (
   select 
      METADATA$FILENAME :: varchar,  
      METADATA$FILE_ROW_NUMBER, 
      $1
   from @vcf_stage 
)
FILE_FORMAT = vcf_stage_fmt 
;

-- Load the Header rows from the Staging table, extracting the Sample_ID from the
-- column header name present in the #CHROM... row of the metadata...
insert into headers (
    filename  ,
    line      ,
    text      
)
Select 
    s.filename,
    s.line,
    s.text
from  vcf_raw_stg s 
where s.text like '#%'
;

-- Load the VCF details from the staging table, again joined to the Sample_ID from the #CHROM... row
-- Materialize a column for each possible attribute in FORMAT field

insert into genotypes (
    CHROM       ,
    POS         ,
    ID          ,
    REF         ,
    ALT         , 
    QUAL        ,
    FILTER      ,
    INFO        ,  
    FORMAT      ,
    GT          ,
    FILENAME    )
with genotype_rows as 
(
  Select 
    substr(split(text, '\t')[0],4)::varchar as chrom, 
    split(text, '\t')[1]::integer as pos, 
    nullif(split(text, '\t')[2],'.')::varchar as id,
    split(text, '\t')[3]::varchar as ref,
    split(split(text, '\t')[4],',')::array as alt,
    nullif(split(text, '\t')[5],'.')::varchar as qual,
    split(text, '\t')[6]::varchar as filter,
    iff(contains(split(text, '\t')[7]::varchar, ';'), split(text, '\t')[7]::varchar, null) as info,
    split(text, '\t')[8]::varchar as format,
    split(split(text, '\t')[9], ':')[array_position('GT'::variant, split(split(text, '\t')[8],':'))]::varchar as gt,
    filename as filename
  from vcf_raw_stg where not text like '#%'
)
select 
    CHROM       ,
    POS         ,
    ID          ,
    REF         ,
    NULLIF(ALT, TO_ARRAY('.')) as ALT    , 
    QUAL        ,
    FILTER      ,
    parse_vcf_info_flat(INFO, get_info_headers_from_stg()) as INFO        ,
    FORMAT      ,
    GT          ,
    FILENAME
from genotype_rows g
order by     
    CHROM ,
    POS   
;

