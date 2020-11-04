-- HELPER UDFs to assist with querying nested elements in VCF files
-- 
use schema GVCF;

-- If using Data Sharing, these UDFs should be created in consumer DBs that have views defined against the share
--   since data sharing does not yet support sharing of UDFs...

-- UDFs to help parse out allele values and reduce query complexity!
CREATE OR REPLACE SECURE FUNCTION allele1 (ref varchar, alt variant, gt varchar)
RETURNS varchar as
$$
case when contains(gt, '|') then 
    case when strtok(gt, '|', 1) = '0' then ref else alt[strtok(gt, '|', 1)::integer-1] end 
else 
    case when strtok(gt, '/', 1) = '0' then ref else alt[strtok(gt, '/', 1)::integer-1] end 
end
$$
;

CREATE OR REPLACE SECURE FUNCTION allele2 (ref varchar, alt variant, gt varchar)
RETURNS varchar as
$$
case when contains(gt, '|') then 
    case when strtok(gt, '|', 2) = '0' then ref else alt[strtok(gt, '|', 2)::integer-1] end 
else 
    case when strtok(gt, '/', 2) = '0' then ref else alt[strtok(gt, '/', 2)::integer-1] end 
end
$$
;

-- Return array search position for a given allele in a stringed pair genotype.  
-- Returns NULL if allele is reference instead of alt.
CREATE OR REPLACE SECURE FUNCTION allele_position (gt varchar, which int)
RETURNS integer as
$$
nullif(
  case when contains(gt, '|') then 
      strtok(gt, '|', which)::integer-1
  else 
      strtok(gt, '/', which)::integer-1
  end,
  -1)
$$
;

-- Boolean function to test if a genotype symbol represents homozygous
CREATE OR REPLACE SECURE FUNCTION is_homozygous (gt varchar)
RETURNS boolean as
$$
  case when contains(gt, '|') then 
      strtok(gt, '|', 1) = strtok(gt, '|', 2)
  else 
      strtok(gt, '/', 1) = strtok(gt, '/', 2) OR strtok(gt, '/', 2) is NULL
  end
$$
;

-- Boolean function to isolate any value in a gvcf at a specific position 
-- E.g. use in the WHERE clause  ... WHERE at_location(pos, end, 23455667) ...
CREATE OR REPLACE SECURE FUNCTION at_location (pos integer, end integer, loc integer)
RETURNS boolean as
$$
    pos <= loc and end >=loc
$$
;

-- Define a join-key for Between Joins linking fixed positions to ranges
CREATE OR REPLACE SECURE FUNCTION genokey (chrom char, pos number) 
RETURNS number as 
$$ 
    (case chrom when 'X' then 23 when 'Y' then 24 else chrom::int end) *1000000000 + pos 
$$ 
;

-- Translate old GT VCF text encoding to one based on a newly assembled ALT array
CREATE OR REPLACE SECURE FUNCTION rebase_gt (old_alt array, new_alt array, gt varchar(5))
RETURNS varchar(5) as
$$
   case when strtok(gt, '/', 1) = '0' 
        then '0'
        else ifnull((array_position (old_alt[strtok(gt, '/', 1)::integer - 1] :: variant, new_alt) +1)  :: varchar, '.') 
   end 
|| 
   iff (contains(gt, '/'), '/' ||
      case when strtok(gt, '/', 2) = '0' 
           then '0'
           else ifnull((array_position (old_alt[strtok(gt, '/', 2)::integer - 1] :: variant, new_alt) +1)  :: varchar, '.')
       end,
       '' ) 
$$
;
