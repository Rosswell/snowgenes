-- A random, evolving collection of sample queries to illustrate use of the 1000Genomes data set...
--
use schema GENOME1K;

-- Alleles --spelled out -- for chr 10, pos 100000-500000, for female samples originating in FINland
select g.sample_id, chrom, pos, ref, alt, gt,
    allele1,
    allele2
from genotypes g
join panel p on p.sample_id = g.sample_id
where chrom = '10' and pos between 100000 and 500000
    and not gt = '0|0' 
    and pop = 'FIN'
    and gender = 'female'
limit 100;


-- Distribution of genotypes in 1kgenome for locations associated with 
-- 'Hereditary_nonpolyposis_colon_cancer' in clinvar CLNDN
-- Note:  Not all of these genotypes shows an allele that is associated with the condition in Clinvar...
select   g.chrom, g.pos, g.ref, g.alt, gt, count (sample_id)
from genotypes g
join clinvar c 
    on c.chrom = g.chrom and c.pos = g.pos and c.ref = g.ref 
where
    array_contains('Hereditary_nonpolyposis_colon_cancer'::variant, c.info:CLNDN)
group by 1,2,3,4,5
order by 1,2,5
;
--> See how the reference counts are ~99% of population

    -- Note that we want to find the samples with actual allele values matching a value in the clinvar ALT 
    -- that indicate a colon cancer association.   Again using the Allele1, Allele2 UDFs to simplify the SQL
    select   g.chrom, g.pos, g.ref, g.alt, gt, count (sample_id)
    from genotypes g
    join clinvar c 
        on c.chrom = g.chrom and c.pos = g.pos and c.ref = g.ref 
        and (
            array_contains(Allele1::variant, c.alt)
            or
            array_contains(Allele2::variant, c.alt)
        )
    where
        array_contains('Hereditary_nonpolyposis_colon_cancer'::variant, c.info:CLNDN)
    group by 1,2,3,4,5
    order by 1,2,5
    ;

-- To zoom-in on variant-specific info, we can flatten the VARIANTS table to the ALLELE level.
-- When flattening, the INDEX is available to extract elements that are allele-specific in the INFO field

-- EXAMPLE:  Variant Alleles with an overall frequency > 60% in Chrom 2
select chrom, pos, ref, alt, v.index alt_allele_id, v.value alt_allele, info:AF[v.index] freq
from variants,
    lateral flatten (INPUT => alt, mode => 'ARRAY') v
    where  info:AF[v.index] > 0.6
    and chrom = '2'
order by 1,2,3,4,5
limit 100
;

-- We get this even more easily using a view VARIANTS_FLATTENED:
select chrom, pos, ref, alt, alt_allele_id, alt_allele, AF freq
from variants_flattened
where af > 0.6
    and chrom = '2'
order by 1,2,3,4,5
limit 100
;

-- LOW-FREQUENCE ALLELE ANALYSIS
-- Show me individuals with <.001 freq allele variants in chrom 11 in a specific 5M bp range 
-- NOTE:  Queries using CTE, Variants_Flattened view, and UDFs: 

        -- Find the individuals will alleles matchhing low-freq variants in the range
        WITH low_freqs as ( 
        select chrom, pos, ref, alt, alt_allele_id, alt_allele, AF freq
        from variants_flattened
            where  AF < 0.001
            and chrom = '11'
            and pos between 65550000 and 70000000
        )
        select v.*, g.sample_id, g.gt, allele1, allele2
        from low_freqs v
            join genotypes g
                on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
            and (v.alt_allele = allele1 OR v.alt_allele = allele2)
        limit 10000;

        -- How many altogether??  Essentially all samples have at least 1 variant in the 5-million bp range!
        WITH low_freqs as ( 
        select chrom, pos, ref, alt, alt_allele_id, alt_allele, AF freq
        from variants_flattened
            where  AF < 0.001
            and chrom = '11'
            and pos between 65550000 and 70000000
        )
        select count(distinct sample_id) from
        (select v.*, g.sample_id, g.gt, allele1, allele2
        from low_freqs v
            join genotypes g
                on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
            and (v.alt_allele = allele1 OR v.alt_allele = allele2)
        ) x
        ;

        -- What's the distribution of variant counts in this range by sample?
        -- How many altogether??
        WITH low_freqs as ( 
        select chrom, pos, ref, alt, alt_allele_id, alt_allele, AF freq
        from variants_flattened
            where  AF < 0.001
            and chrom = '11'
            and pos between 65550000 and 70000000
        )
        select  sample_id, count(*) as cnt from
        (select v.*, g.sample_id, g.gt, allele1, allele2
        from low_freqs v
            join genotypes g
                on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
            and (v.alt_allele = allele1 OR v.alt_allele = allele2)
        ) x
        group by 1 order by 2 desc
        ;        

        -- How many within 4000 bp of locus? -- Notice much smaller set of samples with rare variants
        -- in the narrower range 
        WITH low_freqs as ( 
        select chrom, pos, ref, alt, alt_allele_id, alt_allele, AF freq
        from variants_flattened
            where  AF < 0.001
            and chrom = '11'
            and pos between 68786233 and 68790000
        )
        select count(distinct sample_id) from
        (select v.*, g.sample_id, g.gt, allele1, allele2
        from low_freqs v
            join genotypes g
                on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
            and (v.alt_allele = allele1 OR v.alt_allele = allele2)
        ) x
        ;

        -- And the narrow list:
        WITH low_freqs as ( 
        select chrom, pos, ref, alt, alt_allele_id, alt_allele, AF freq
        from variants_flattened
            where  AF < 0.001
            and chrom = '11'
            and pos between 68786233 and 68790000
        )
        select v.*, g.sample_id, g.gt, allele1, allele2
        from low_freqs v
            join genotypes g
                on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
            and (v.alt_allele = allele1 OR v.alt_allele = allele2)
            order by chrom, pos
        ;

-- Use thge Materialized View (genotypes-by-sample) ordered by Sample to find all the rare homozygous sites in an individual:
-- Find **all** the low freq homozygous alleles present in sample HG02215
        WITH low_freqs as ( 
        select chrom, pos, ref, alt, alt_allele_id, alt_allele, AF freq
        from variants_flattened
            where  AF < 0.001
        )
        select * from low_freqs v
            join genotypes_by_sample g
                on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
                and (v.alt_allele = allele1 OR v.alt_allele = allele2)
            where sample_id = 'HG02215' and is_homozygous(gt)
            order by g.chrom, g.pos
        ;


-- Big Query Range Join, Allele Frequency example 
-- from https://cloud.google.com/genomics/docs/how-tos/interval-joins 
--  
-- (Snowflake <20 secs against Genomes1K, vs. BQ taking 90 secs!  BQ charges for 3 TB scan ($15) vs. effective Snow Pruning!))
-- Note:  we use SNP annotation DB instead of the TUTE table which we don't have.  Identically challenging :
with intervals as (
  select 'PRCC' as gene, '1' as chr, 156736274 as gene_start, 156771607 as gene_end, 156636274 as region_start, 156871607 as region_end from dual
  UNION
  select 'NTRK1' as gene, '1' as chr, 156785541 as gene_start, 156852640 as gene_end, 156685541 as region_start, 156952640 as region_end from dual
  UNION
  select 'PAX8' as  gene, '2' as chr, 113972574 as gene_start, 114037496 as gene_end, 113872574 as region_start, 114137496 as region_end from dual
  UNION
  select 'FHIT' as  gene, '3' as chr, 59734036 as gene_start, 61238131 as gene_end, 59634036 as region_start, 61338131 as region_end from dual
  UNION
  select 'PPARG' as gene, '3' as chr, 12328349 as gene_start, 12476853 as gene_end, 12228349 as region_start, 12576853 as region_end from dual
),
my_flattened_variants as (
    select 
        v.CHROM, v.POS, v.ref, v.alt, v.alt_allele, i.gene, v.AF,
        sum(iff(v.alt_allele = allele1, 1, 0) + iff(v.alt_allele = allele2, 1, 0)) as num_variant_alleles,
        sum(allele_count(gt)) as total_num_alleles
    from variants_flattened v
    join genotypes g
         on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
    join intervals i
        on v.chrom = i.chr and i.region_start <= v.pos and i.region_end >= v.pos
  where v.AF < 0.01
  group by  v.CHROM, v.POS, v.ref, v.alt, v.alt_allele, i.gene, v.AF
  )
select v.*, c.alt, c.info from my_flattened_variants v
left join snp c on v.chrom = c.chrom and v.pos = c.pos and v.ref = c.ref 
and array_contains(v.alt_allele, c.alt)
order by 1,2,5
-- limit 100
;  

--  CPIC® Guideline for Clopidogrel and CYP2C19  
--  https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/
--  Use SNP annotation DB to find rare variants for a specific gene :
with my_flattened_variants as (
    select 
        v.CHROM, v.POS, v.ref, v.alt, v.alt_allele, v.AF, s.info:GENEINFO geneinfo, v.id,
        sum(iff(v.alt_allele = allele1, 1, 0) + iff(v.alt_allele = allele2, 1, 0)) as num_variant_alleles,
        sum(allele_count(gt)) as total_num_alleles
    from variants_flattened v
    join genotypes g
         on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
    join snp s
        on v.chrom = s.chrom and s.pos = v.pos and s.alt = v.alt
  where s.info:GENEINFO like  'CYP2C19%' 
  and v.AF < 0.01
  group by  1,2,3,4,5,6,7,8
  )
select * from my_flattened_variants v
;

-- Or drill down using specific RS numbers in the literature 
-- e.g. https://api.pharmgkb.org/v1/download/file/attachment/CYP2C19_allele_definition_table.xlsx
with my_flattened_variants as (
    select 
        v.CHROM, v.POS, v.ref, v.alt, v.alt_allele, v.AF, s.info:GENEINFO geneinfo, v.id,
        sum(iff(v.alt_allele = allele1, 1, 0) + iff(v.alt_allele = allele2, 1, 0)) as num_variant_alleles,
        sum(allele_count(gt)) as total_num_alleles
    from variants_flattened v
    join genotypes g
         on g.chrom = v.chrom and g.pos = v.pos and g.ref = v.ref and g.alt = v.alt
    join snp s
        on v.chrom = s.chrom and s.pos = v.pos and s.alt = v.alt
  where s.info:GENEINFO like  'CYP2C19%' 
  and v.AF < 0.01
  group by  1,2,3,4,5,6,7,8
  )
select * from my_flattened_variants v
where v.id[0] in (
'rs11188072',
'rs113164681',
'rs111490789',
'rs17878739',
'rs7902257',
'rs12248560',
'rs367543001',
'rs28399504',
'rs367543002',
'rs367543003',
'rs55752064',
'rs17882687',
'rs17885098',
'rs145328984',
'rs118203756',
'rs12769205',
'rs41291556',
'rs72552267 ',
'rs17884712',
'rs58973490',
'rs140278421',
'rs370803989',
'rs4986893',
'rs6413438',
'rs4244285',
'rs72558186',
'rs138142612',
'rs3758580',
'rs3758581',
'rs118203757',
'rs113934938',
'rs377184510',
'rs17879685',
'rs17886522',
'rs56337013',
'rs192154563',
'rs118203759',
'rs55640102') 
;
