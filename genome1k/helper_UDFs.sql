-- HELPER UDFs to assist with querying nested elements in VCF files
-- 
use schema GENOME1K;

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

-- Boolean function to determine whether Allele1 or Allele2 is contained in a target array
CREATE OR REPLACE SECURE FUNCTION allele_contained_in (ref varchar, alt variant, gt varchar, target variant)
RETURNS boolean as
$$
    array_contains(Allele1(ref, alt, gt)::variant, target)
    or
    array_contains(Allele2(ref, alt, gt)::variant, target)
$$
;

-- Boolean function to determine whether Allele1 or Allele2 is equal to a target value
CREATE OR REPLACE SECURE FUNCTION allele_matches (ref varchar, alt variant, gt varchar, target varchar)
RETURNS boolean as
$$
    target = allele1(ref, alt, gt) 
    or 
    target = allele2(ref, alt, gt)
$$
;

-- The count of alleles reflected in a symbolic GT value
CREATE OR REPLACE SECURE FUNCTION allele_count (gt varchar)
RETURNS integer as
$$
    regexp_count(gt, '[/|]') + 1
$$
;

-- Number of alleles at a position that match a target value 
-- Useful primitive for creating distribution statistics across a population
CREATE OR REPLACE SECURE FUNCTION matching_allele_count (allele_to_match varchar, ref varchar, alt variant, gt varchar)
RETURNS integer as
$$
    iff(allele_to_match = allele1(ref, alt, gt),1,0) + 
    iff(allele_to_match = allele2(ref, alt, gt),1,0)
$$
;