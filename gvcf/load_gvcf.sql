-- ALTERNATIVE Loading script for gVCF files
-- Assumes that gVCF is pre-processed to generate pairs of files reflecting headers and detail rows:
--    *.vcfrows.gz
--    *.vcfheaders.gz
-- both in text format with gz compression
-- 
-- This loading approach is designed for faster loads but requires some preprocessing
--
-- The final detail VCF data, rather than having separate columns materialized for every FORMAT field,
--    instead relies on a VIEW definition to parse-out less frequently used fields on demand, reducing load
--    overhead.
-- 
-- The load process is as follows:
--   Headers are loaded into staging and then moved into permanent table with a SAMPLE_ID lookup
--   Detail is loaded directly as an INSERT from staging joined to SAMPLE_ID lokup based on header 
--
use schema GVCF;

CREATE OR REPLACE STAGE gvcfstage
    URL = 's3://****/'
    CREDENTIALS = (  AWS_KEY_ID =  '*****' AWS_SECRET_KEY = '*****' )
    FILE_FORMAT = (TYPE = CSV FIELD_DELIMITER = '\t' COMPRESSION = gzip ERROR_ON_COLUMN_COUNT_MISMATCH = FALSE);

-- format for vcf headers
 create or replace file format gvcf_header_fmt  TYPE = CSV COMPRESSION = AUTO FIELD_DELIMITER = '`';

copy into HEADERS (
    FILENAME  ,
    SAMPLE_ID ,
    line      ,
    text      
)  from (
   select 
      split(METADATA$FILENAME,'/')[1] :: varchar,  
      split(split(METADATA$FILENAME,'/')[1],'.')[0] :: varchar,
      METADATA$FILE_ROW_NUMBER, 
      $1
   from @gvcfstage/gvcf/
)
FILE_FORMAT = gvcf_header_fmt 
PATTERN =   '.*.vcfheaders.gz'
;

copy into genotypes_raw (
    CHROM       ,
    POS         ,
    END         ,
    ID          ,
    REF         ,
    ALT_STRING  , 
    QUAL        ,
    FILTER      ,
    INFO        ,
    SAMPLE_ID   ,
    FORMAT      ,
    GT          ,
    VALS        ,
    FILENAME    
) from (
    Select
        substr($1,4) as CHROM,
        $2 as POS,
        ifnull(split($8, 'END=')[1], $2) as END,
        nullif($3, '.') as ID,
        $4 as REF,
        $5 as ALT_STRING,
        nullif($6, '.') as QUAL,
        $7 as FILTER,
        iff(startswith($8,'END='), '.', $8) as INFO,
        strtok(strtok(METADATA$FILENAME, '/', 2), '.', 1) as SAMPLE_ID,
        $9 as FORMAT,
        strtok($10, ':', 1) as GT,
        $10  as VALS,
        strtok(METADATA$FILENAME, '/', 2) as FILENAME 
  from @gvcfstage/gvcf/
)
PATTERN = '.*.vcfrows.gz'
ON_ERROR = 'skip_file'
;

-- load clinvar data for hg38 format
copy into clinvar_raw from @gvcfstage
    files = ('clinvar.vcf.gz')
    FILE_FORMAT = (TYPE = CSV FIELD_DELIMITER = '\t' COMPRESSION = auto SKIP_HEADER = 28)
    ON_ERROR = SKIP_FILE
;

-- load clinvar vcf header
insert into headers (filename, sample_id, line, text) 
    select 
    'clinvar',
    NULL,
    METADATA$FILE_ROW_NUMBER, 
    $1
    from @gvcfstage/clinvar.vcf.gz   (FILE_FORMAT => gvcf_header_fmt ) t
    where $1 like '#%'
;

insert into CLINVAR (
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
    split(id, ';') as id, 
    ref, 
    split(alt, ',') as alt, 
    qual, 
    filter, 
    parse_vcf_info(info, get_info_headers('clinvar')) INFO
from CLINVAR_RAW
order by 1,2,4,5
;


-- load snp data for hg38 format
copy into snp_raw from @gvcfstage 
    PATTERN = '.*snp/part.*'
    ON_ERROR = SKIP_FILE
;

-- load snp vcf header
copy into headers (filename, sample_id, line, text) from (                                
        select 
        'snp',
        NULL,
        METADATA$FILE_ROW_NUMBER, 
        $1
    from @gvcfstage
) FILE_FORMAT = ( FORMAT_NAME = 'gvcf_header_fmt' )
  PATTERN = '.*snp/header.*'
; 

insert into SNP (
  CHROM,
  ORIG_CHROM,
  POS,
  ID,
  REF,
  ALT,
  QUAL,
  FILTER,
  INFO
)
select
    coalesce(c.chrom, g.chrom) CHROM,
    g.chrom ORIG_CHROM,
    pos, 
    split(id, ';') as id, 
    ref, 
    split(alt, ',') as alt, 
    qual, 
    filter, 
    parse_vcf_info(info, get_info_headers('snp')) INFO
from SNP_RAW g
left outer join GENBANK_TO_CHROM c
    on g.CHROM = c.GENBANK_ID
order by 1,2,4,5
;

-- format for BED files
 create or replace file format gvcf_bed_fmt  TYPE = CSV COMPRESSION = AUTO FIELD_DELIMITER = '\t';


copy into EXAC_BED 
from @gvcfstage/bed/
)
FILE_FORMAT = gvcf_bed_fmt
files =   'exac.bed'
;

create or replace table exac_bed_positions as
    select chrom, start_pos + offset as pos 
    from exac_bed cross join offsets 
    where end_plus1 > start_pos + offset
    order by chrom, pos
;
