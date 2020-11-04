-- SQL Script to load data into GNOMAD schema
use schema GNOMAD;

truncate table variants_raw;
truncate table variants;
truncate table headers;

-- Stage to hold both headers and tab-delimited VCF body lines
create or replace stage gnomadstage 
    URL = 's3://sfc-benchmarks/gnomad/'
    CREDENTIALS = (  AWS_KEY_ID =  '*************' AWS_SECRET_KEY = '************' )
    FILE_FORMAT = (TYPE = CSV FIELD_DELIMITER = '\t' COMPRESSION = gzip);

-- format for vcf headers
 create or replace file format header_fmt  TYPE = CSV COMPRESSION = AUTO FIELD_DELIMITER = '`';

-- Load the chromosome header files for EXOME
copy into headers (vcfname, line, text) from (                                
        select 
        regexp_replace (METADATA$FILENAME,
        '.*gnomad\\/(.*)\\/.*',
        '\\1'), 
        METADATA$FILE_ROW_NUMBER, 
        $1
    from @gnomadstage
) FILE_FORMAT = ( FORMAT_NAME = 'header_fmt' )
  PATTERN = '.*exome.*chr.*/header.*'
; 

-- Load the raw VCF data for EXOME.  
-- Load files from chr1-chr22, X and Y
copy into variants_EXOME_raw from @gnomadstage 
    PATTERN = '.*exome.*chr.*/part.*'
    ON_ERROR = SKIP_FILE
;

insert into VARIANTS_EXOME (
  CHROM,
  POS,
  ID,
  REF,
  ALT,
  QUAL,
  FILTER,
  INFO
)
select 
    chrom, 
    pos, 
    id, 
    ref, 
    alt, 
    qual, 
    filter, 
    parse_vcf_info_flat(info, get_info_headers()) INFO
from VARIANTS_EXOME_RAW
order by 1,2,4,5
;

-- Load the chromosome header files for GENOME
copy into headers (vcfname, line, text) from (                                
        select 
        regexp_replace (METADATA$FILENAME,
        '.*gnomad\\/(.*)\\/.*',
        '\\1'), 
        METADATA$FILE_ROW_NUMBER, 
        $1
    from @gnomadstage
) FILE_FORMAT = ( FORMAT_NAME = 'header_fmt' )
  PATTERN = '.*genome.*chr.*/header.*'
; 

-- Load the raw VCF data for GENOME.  
-- Load files from chr1-chr22, X and Y
copy into variants_GENOME_raw from @gnomadstage 
    PATTERN = '.*genome.*chr.*/part.*'
    ON_ERROR = SKIP_FILE
;

insert into VARIANTS_GENOME (
  CHROM,
  POS,
  ID,
  REF,
  ALT,
  QUAL,
  FILTER,
  INFO
)
select 
    chrom, 
    pos, 
    id, 
    ref, 
    alt, 
    qual, 
    filter, 
    parse_vcf_info_flat(info, get_info_headers()) INFO
from VARIANTS_GENOME_RAW
order by 1,2,4,5
;
