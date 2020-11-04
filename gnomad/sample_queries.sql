-- Example Queries using GNOMAD summary data joined to 1KGenomes 
use schema GNOMAD;

-- Identify variants of interest based on statistical properties in GNOMAD
--    E.g. Finland-origin frequency > 0.5 but overall pop freq < 0.1
select info:AF_fin, info:AF, * 
from variants_exome
where info:AF_fin > 0.5 and info:AF < 0.1;

-- Compare those positions with overall stats in 1kGenomes variants
with targets as (
  select info:AF_fin, info:AF, 
  chrom, pos, alt_allele, ref 
  from variants_exome 
  where info:AF_fin > 0.5 and info:AF < 0.1)
select t.*, v.* from genome1k.variants_flattened v
join targets t
    on v.chrom = t.chrom
    and v.pos = t.pos
    and t.alt_allele = v.alt_allele
order by t.chrom, t.pos
;

-- Find counts of individuals in 1kGenomes that match
with targets as (
  select info:AF_fin, info:AF, 
  chrom, pos, alt_allele, ref 
  from variants_exome 
  where info:AF_fin > 0.5 and info:AF < 0.1)
select 
    t.*, count(*) as cnt
from 
    genome1k.genotypes g
join targets t
on  g.chrom = t.chrom
    and g.pos = t.pos
    and (t.alt_allele = genome1k.allele1(g.ref, g.alt, g.gt) OR t.alt_allele = genome1k.allele2(g.ref, g.alt, g.gt))
group by 1,2,3,4,5,6
order by chrom, pos;

-- Now let's compare to the country of origin from the 1kGenomes pop data.  Are we predictive???
with targets as (
  select info:AF_fin, info:AF, 
  chrom, pos, alt_allele, ref 
  from variants_exome
  where info:AF_fin > 0.5 and info:AF < 0.1)
select 
    t.*, p.pop, count(*) as cnt
from 
    genome1k.genotypes g
join targets t
on  g.chrom = t.chrom
    and g.pos = t.pos
    and (t.alt_allele = genome1k.allele1(g.ref, g.alt, g.gt) OR t.alt_allele = genome1k.allele2(g.ref, g.alt, g.gt))
join genome1k.panel p
    on p.sample_id = g.sample_id 
group by 1,2,3,4,5,6,7
order by cnt desc;