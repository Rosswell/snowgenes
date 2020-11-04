-- EXAMPLES:
use schema GVCF;

-- Get all the genotypes AT a position
select * from genotypes_by_position -- chr 2 pos 699000
where 
    chrom = '2' and 
    pos <= 700000 and end >=700000
; 

select * from genotypes_by_position -- chr11 pos 64497189
where 
    chrom = '11' and 
    pos <= 64497189 and end >=64497189
; 

-- Note above that ref-values for range entries are not known until we find a variant defined in the data...


-- Non-refs (variants) in a range that PASS a filter:
select pos, end, any_value(ref), count(sample_id) 
from genotypes_by_position
where 
  chrom = '11' 
  and pos between 64000000 and 64500000
  and filter='PASS'
  and not GT in ('0/0', './.', '0', '.')
group by pos, end
order by pos
limit 1000;

-- Get all samples' genotype values for observed variant positions in the range
with nonrefs as (
    select chrom, pos, end, any_value(ref) as ref
    from genotypes_by_position
    where 
        chrom = '11' 
        and pos between 64000000 and 64500000
        and filter='PASS'
        and not GT in ('0/0', './.', '0', '.')
    group by chrom, pos, end
)
select 
    sample_id,
    n.chrom,
    n.pos as pos,
    n.ref,
    filter,
    g.pos gpos,
    g.end gend,
    g.gt,
    case when strtok(gt, '/', 1) = '0' then n.ref else alt[strtok(gt, '/', 1)::integer-1] end as allele1,
    case when strtok(gt, '/', 2) = '0' then n.ref else alt[strtok(gt, '/', 2)::integer-1] end as allele2,
    dp,
    gq,
    min_dp
from genotypes_by_position g
join nonrefs n
    on n.chrom = g.chrom and n.pos >= g.pos and n.end <= g.end
-- where filter='PASS'
order by 2,3,1
--limit 1000
;


-- What are the counts of above genotype values per position??
with nonrefs as (
    select chrom, pos, end, any_value(ref) as ref
    from genotypes_by_position
    where 
        chrom = '11' 
        and pos between 64000000 and 64500000
        and filter='PASS'
        and not GT in  ('0/0', './.', '0', '.')
    group by chrom, pos, end
)
select chrom, pos, ref, allele1, allele2, gt, count(*) from (
select 
    sample_id,
    n.chrom,
    n.pos as pos,
    n.ref,
    filter,
    g.pos gpos,
    g.end gend,
    g.gt,
    case when split(gt,'/')[0] = '0' then n.ref else alt[split(gt,'/')[0]::integer-1] end as allele1,
    case when split(gt,'/')[1] = '0' then n.ref else alt[split(gt,'/')[1]::integer-1] end as allele2,
    dp,
    gq,
    min_dp
from genotypes_by_position g
join nonrefs n
    on n.chrom = g.chrom and n.pos >= g.pos and n.end <= g.end) x
group by 1, 2, 3, 4, 5, 6
order by 1, 2
;

-- Construct an ALT array as the UNION of all ALT values in a range, per position
-- As a first step, observe all the alt values per position
select chrom, pos, end, ref, f.value
from genotypes_by_position g,
    lateral flatten(input => g.alt) f
where 
    chrom = '11' 
    and pos between 64000000 and 64500000
    and filter='PASS'
    and not GT in  ('0/0', './.', '0', '.')
;

-- Here is how we Dynamically create an ALT array for a position based on all ALTs in rows for that position
select distinct chrom, pos, end, ref, array_agg(distinct f.value) as alt
from genotypes_by_position g,
    lateral flatten(input => g.alt) f
where 
    chrom = '11' 
    and pos between 64000000 and 64500000
    and filter='PASS'
    and not GT in ('0/0', './.', '0', '.')
group by 1,2,3,4
;

-- We use this to package a VCF with freshly-defined ALT values -- usable even when samples at a position
-- utilize different values or ordering of ALTs
with alts as (
    select distinct chrom, pos, end, ref, any_value(id) as id, array_agg(distinct f.value) as alt
    from genotypes_by_position g,
        lateral flatten(input => g.alt) f
    where 
        chrom = '11' 
        and pos between 64000000 and 64500000
        and filter='PASS'
        and not GT in ('0/0', './.', '0', '.')
    group by 1,2,3,4 
),
results as (
select 
    sample_id,
    n.chrom as CHROM,
    n.pos as POS,
    n.id as ID,
    n.ref as REF,
    ARRAY_TO_STRING(n.alt, ',') as ALT,
    '.' as QUAL,
    '.' as FILTER,
    '.' as INFO,
    'GT:DP:GQ:MIN_DP:FT' as FORMAT,
    case when strtok(gt, '/', 1) = '0' then 0 :: varchar else ifnull((array_position (g.alt[strtok(gt, '/', 1)::integer - 1] :: variant, n.alt) +1)  :: varchar, '.') end || 
      iff (contains(gt, '/'), 
           '/' || case when strtok(gt, '/', 2) = '0' then 0 :: varchar else ifnull((array_position (g.alt[strtok(gt, '/', 2)::integer - 1] :: variant, n.alt) +1)  :: varchar, '.') end,
           '' ) 
    || ':' || coalesce(dp :: varchar, '.')
    || ':' || coalesce(gq :: varchar, '.')
    || ':' || coalesce(min_dp :: varchar,'.') 
    || ':' || coalesce(g.filter, '.')  as val
from genotypes_by_position g
join alts n
    on n.chrom = g.chrom and n.pos >= g.pos and n.end <= g.end
)
select * from results
pivot(max (val)
     for sample_id in ('SlPFTOJi', 'RcDfJotP', 'XTsSvboH', 'JeBTaIMX', 'fAeYoHXC', 'UgyTTpve', 'kjicuIpB', 'mCGYrLln', 'ZWZLUyfm', 'VHMpiTYC', 'mKHvZFAJ', 'bXDtVeuv', 'xzrOKjYV', 'FdUwKlXj', 
                       'NFiQejWN', 'cLHWhvEw', 'tEihALIR', 'hFmIuoPN', 'wkFmnQbK', 'KfciHznG', 'weTnPgoq', 'EmGpcBww', 'cUuAWNsR', 'SGnPTIBV', 'UKRJPuwg', 'zkQHvQhJ', 'xRAODqsm', 'lHrvGfVv', 
                       'rFTYmDps', 'BZxiLVPf', 'rduTsxva', 'pFtYwIBX', 'IIdaEdZV', 'tIySMRDw', 'NZChVRXi', 'pffFKYvX', 'pLEKZEIr', 'dUwKFaRn', 'PiplAHTj', 'QCndDLsD', 'dBbPdUlA', 'wqHIJbDa', 
                       'TNbDlHCf', 'zfdgjsRb', 'ETcWRDHM', 'lMwBGhDw', 'ACHUcbBa', 'uJXyOXlc', 'qCMVpZQl', 'SlYSaFui', 'DJxhiGZC', 'XNDkAIqZ', 'sFrDtfQb', 'NlsJlrzR', 'XnujMpHG', 'bOSljHhM', 
                       'YbBcsxtA', 'ZtOfNcje', 'zuEpXUol', 'GIZiSPPw', 'cZYKwIDU', 'reIsriqy', 'HspNyazF', 'xIomjyOk', 'QdtLdnrp', 'bHnFyTXn', 'cGAJnKtS', 'LiDDnpya', 'TRjMHaFJ', 'RXtghUJC', 
                       'RtdjvUrW', 'HJPRfZgG', 'hXXFpiwY', 'soDozTcS', 'DPuMoNgV', 'FtzJARCt', 'aaFvQRFt', 'kGoUEvXN', 'FsRmytGF', 'FgdHPLgA', 'zughZyzZ', 'iVESyEBj', 'uBUqVVXL', 'IVHBijAE', 
                       'bvmGVeNH', 'PLkSbyZE', 'NhsOoajx', 'eaHJdfvY', 'jKWbZNLQ', 'EGBjwfcQ', 'WuIlJmzh', 'vDyKgiIS', 'ZoJfcIuK', 'nkBnbEyZ', 'zCxhDrQS', 'ZZduePKn', 'AHszUkGn', 'wmNohzLL', 
                       'yCwQufVW', 'wOoSvvXx')) 
     as p 
order by 1,2
;



-- Data Validation Queries
with alt_variations as (
select chrom, pos, ref, count(distinct alt) as varcnt
from genotypes_by_position
  where pos = end and filter='PASS' and not GT in ('0/0', './.', '0', '.')
group by chrom, pos, ref)
select * from alt_variations where varcnt > 1
limit 1000;

with alt_refs as (
select chrom, pos, count(distinct ref) as refcnt
from genotypes_by_position
  where pos = end and filter='PASS'
group by chrom, pos)
select * from alt_refs where refcnt > 1
limit 1000;

-- Occasional Anomalous REF values:  how to handle?
select * from 
genotypes_by_position
where chrom = '6' and pos = 30888161 and pos = end
;

select * from 
genotypes_by_position
where chrom = '11' and pos = 1018562 and pos = end
;

-- Queries to select data subset to use in analyzing data-defined variants
-- First create staging tables for the analysis:

-- Variants derived from cohort:
select count(*) as c, chrom, pos, ref, alt_STRING, gt
from GENOTYPES_RAW
where ALT_STRING != '<NON_REF>'
--and filter = 'PASS'
--and sample_id between  'v1_C' and 'v1_D'
group by 2,3,4,5,6
order by 2,3;

-- Depth distribution:
select 
    count(*) as c, 
    (dp/5)::integer * 5 as depth
from GENOTYPES_BY_SAMPLE
where 
ALT_STRING != '<NON_REF>'
    and filter = 'PASS'
group by 2
order by 2;

-- Filter distribution:
select 
   count(*) as count, 
   filter
from GENOTYPES_BY_SAMPLE
group by 2
order by 2;


-- Full cohort frequencies for a subset of the genome using a BED file
-- Chrom 1 pos 1-1000000 - 3 secs
with myGenotypes as (
   select  chrom, pos, end, ref, alt_string, gt
   from genotypes_by_sample
   where chrom = '1' and pos between 1 and 1000000),
myPositions as (
   select chrom, pos 
   from exac_bed_positions
   where chrom = '1' and pos between 1 and 1000000)
select count(*) as c, b.chrom, b.pos, g.alt_string, g.gt
from myGenotypes g, myPositions b
   where g.pos <= b.pos
   and g.end >= b.pos
group by 2,3,4,5;

--*******************************************************************************************************
--  CPIC® Guideline for Clopidogrel and CYP2C19  
--  https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/
-- Use HG38 coordinates of variants of interest to find the multi-allelic variants that are candidates:

CREATE or replace TABLE CYP2C19_alleles (chrom varchar(2), pos int, ref char(1), alt char(1));

INSERT INTO CYP2C19_alleles(Chrom, Pos, ref, alt) values
    ('10', 94759304, 'C', 'T'),
    ('10', 94760676, 'C', 'T'),
    ('10', 94760686, 'C', 'A'),
    ('10', 94761267, 'T', 'C'),
    ('10', 94761665, 'G', 'A'),
    ('10', 94761900, 'C', 'T'),
    ('10', 94762693, 'G', 'A'),
    ('10', 94762706, 'A', 'G'),
    ('10', 94762712, 'C', 'T'),
    ('10', 94762715, 'T', 'C'),
    ('10', 94762755, 'T', 'C'),
    ('10', 94762760, 'A', 'C'),
    ('10', 94762788, 'A', 'T'),
    ('10', 94762804, 'C', 'T'),
    ('10', 94762856, 'A', 'G'),
    ('10', 94775106, 'C', 'T'),
    ('10', 94775121, 'C', 'T'),
    ('10', 94775160, 'G', 'C'),
    ('10', 94775185, 'A', 'G'),
    ('10', 94775367, 'A', 'G'),
    ('10', 94775416, 'T', 'C'),
    ('10', 94775453, 'G', 'A'),
    ('10', 94775489, 'G', 'A'),
    ('10', 94775507, 'G', 'A'),
    ('10', 94780574, 'G', 'C'),
    ('10', 94780579, 'G', 'A'),
    ('10', 94780653, 'G', 'A'),
    ('10', 94781858, 'C', 'T'),
    ('10', 94781859, 'G', 'A'),
    ('10', 94781944, 'G', 'A'),
    ('10', 94781999, 'T', 'A'),
    ('10', 94842861, 'G', 'A'),
    ('10', 94842865, 'C', 'T'),
    ('10', 94842866, 'A', 'G'),
    ('10', 94842879, 'G', 'A'),
    ('10', 94842995, 'G', 'A'),
    ('10', 94849964, 'A', 'G'),
    ('10', 94849995, 'C', 'T'),
    ('10', 94850018, 'A', 'C'),
    ('10', 94852738, 'C', 'T'),
    ('10', 94852765, 'C', 'T'),
    ('10', 94852785, 'C', 'G'),
    ('10', 94852914, 'A', 'C')
;

-- Unfortunately in Syntheticv gvcf the altcounts are all 0, but this query only takes 30 secs against the 
-- MV organized by chrom, pos -- doing all the heavy lifting anyway
select multi_genotype, altcount, count(*) cnt from (
  Select sample_id, 
  listagg(coalesce(allele1(a.ref, g.alt, g.gt), '.')||'-'||coalesce(allele2(a.ref, g.alt, g.gt), '.'), ' , ') WITHIN GROUP (order by a.pos) multi_genotype,
  sum(case when (a.ref != allele1(a.ref, g.alt, g.gt) and allele1(a.ref, g.alt, g.gt) IS NOT NULL) 
      OR (a.ref != allele2(a.ref, g.alt, g.gt) and allele2(a.ref, g.alt, g.gt) IS NOT NULL) then 1 else 0 end) altcount
  from genotypes_by_position g join CYP2C19_alleles a
       on a.pos >= g.pos and a.pos <=g.end
  where a.chrom = '10' and g.chrom = '10'
  and g.pos <= 94852914 and g.end >= 94759304 -- narrow range further
  group by sample_id) x
group by 1, 2
order by 2,3 desc
;

--*****************************************************************************
-- Process for generating multisample VCF -- using a small range as example:
-- Longest step: 2 mins on 4XL --> this filtered base table can remain semi-permanent
create or replace transient table Genotypes_1_100000 as
select * from GENOTYPES_BY_SAMPLE
where chrom = '1' and pos < 1000000;

-- Remaining steps are a few secs each on 4XL
-- Variant positions and all genotypes at them identified
create or replace temporary table my_Variants_1_1000000 as 
	select chrom, pos, max(id) as id, ref, alt, gt, count(*) as cnt
	from Genotypes_1_100000
	where ALT_STRING != '<NON_REF>'
	group by 1,2,4,5,6
    order by 1 asc,2 asc,7 desc;

-- Redefine the ALT field to be a common list for a given position
create or replace temporary table MY_VARIANT_LOCI_1_1000000 as 
    select chrom, pos, id, ref, alt, array_to_string(alt, ',') as alt_string, cnt 
    from (
      select chrom, pos, id, ref, array_agg(distinct f.value) as alt, sum(cnt) as cnt
      from MY_VARIANTS_1_1000000 v,
          lateral flatten(input => v.alt) f
      group by 1,2,3,4
      order by 1,2
    )x
;

-- detail results for VCF - LIST
create or replace temporary table export_detail_calls_list_1_1000000 as
        with myGenotypes as (
        select  chrom, pos, end, alt, gt, sample_id
        from Genotypes_1_100000
        where chrom = '1'
        ),
        myPositions as (
        select chrom, pos, ref, alt, alt_string 
        from MY_VARIANT_LOCI_1_1000000
        where chrom = '1'
        )
        select  b.chrom, b.pos, b.ref, b.alt_string, 
        rebase_gt(g.alt, b.alt, g.gt) as gt, 
        listagg(sample_id,',')  WITHIN GROUP (order by sample_id) as samples,
        count(*) as cnt
        from myGenotypes g, myPositions b
        where g.pos <= b.pos
        and g.end >= b.pos
        group by 1,2,3,4,5
        order by 1,2;

select * from export_detail_calls_list_1_1000000 order by 1,2 limit 100;

****************************************************************************