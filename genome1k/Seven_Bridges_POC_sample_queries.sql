-- QUERY 1 - BigQuery:
------------------------
WITH filtered_snp_calls AS (
     SELECT  reference_name, c.call_set_name name, CONCAT(reference_bases, '->', alternate_bases[ORDINAL(1)]) AS mutation
     FROM  `genomics-public-data.1000_genomes_phase_3.variants` AS v,
                  UNNEST(v.call) AS c
     WHERE
     # Only include biallelic SNPs.
     reference_bases IN ('A','C','G','T')
     AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
     AND (ARRAY_LENGTH(alternate_bases) = 1
     OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)] = '<*>'))
     # Skip homozygous reference calls and no-calls.
     AND EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g > 0)
     AND NOT EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g < 0)
     ),
     mutation_type_counts AS (
     SELECT
     reference_name,
     name,
     SUM(CAST(mutation IN ('A->G', 'G->A', 'C->T', 'T->C') AS INT64)) AS transitions,
     SUM(CAST(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
     'A->T', 'T->A', 'C->G', 'G->C') AS INT64)) AS transversions
     FROM filtered_snp_calls
     GROUP BY
     reference_name,
     name
     )
     SELECT reference_name,  name, transitions, transversions, transitions/transversions AS titv
     FROM mutation_type_counts
     WHERE transversions > 0
     ORDER BY titv DESC,  name;


-- QUERY 2 - BigQuery:
------------------------
with gcount as (
  SELECT v.start, v.reference_bases, TO_JSON_STRING(v.alternate_bases) alt_bases, 
           sum(if( c.genotype[OFFSET(0)] != -1 and (c.genotype[OFFSET(0)] != 0 or c.genotype[OFFSET(1)] != 0),1,0)) called_non_ref,
           sum(if( c.genotype[OFFSET(0)] != -1 and (c.genotype[OFFSET(0)] != c.genotype[OFFSET(1)]),1,0)) called_het,
           sum(if( c.genotype[OFFSET(0)] != -1,1,0)) called_cnt,
           sum(if( c.genotype[OFFSET(0)] = -1,1,0)) not_called_cnt 
  FROM `genomics-public-data.1000_genomes_phase_3.variants` as v,
  UNNEST(v.call) AS c
  where v.reference_name = '20'
  group by v.start, v.reference_bases, alt_bases
),
homozyg_count as (
  with prep as (
     select v.start, v.reference_bases, TO_JSON_STRING(v.alternate_bases) alt_bases, TO_JSON_STRING(c.genotype) gtype, count(TO_JSON_STRING(c.genotype)) counter
         FROM `genomics-public-data.1000_genomes_phase_3.variants` as v,
               UNNEST(v.call) AS c
         where v.reference_name = '20'
               and c.genotype[OFFSET(0)] = c.genotype[OFFSET(1)]
         group by v.start, v.reference_bases, alt_bases, gtype
         order by v.start, v.reference_bases, alt_bases, gtype
 )
  select g.start, g.reference_bases, g.alt_bases, TO_JSON_STRING(ARRAY_AGG(STRUCT(g.gtype as genotype, counter as count))) homozyg
     from prep as g
     group by g.start, g.reference_bases, g.alt_bases
     order by g.start, g.reference_bases, g.alt_bases
)
select v.start, v.end, v.reference_bases, v.alternate_bases,
       STRUCT(h.homozyg as homozygote_count,
                        g.called_cnt/v.AN*2.0 as call_rate,
                        g.called_cnt as n_called,
                        g.not_called_cnt as n_not_called,
                        g.called_het as n_het,
                        g.called_non_ref as n_non_ref) as variant_qc
    from `genomics-public-data.1000_genomes_phase_3.variants` as v,
         gcount g,
         homozyg_count h
    where g.start = v.start
          and g.reference_bases = v.reference_bases
          and g.alt_bases = TO_JSON_STRING(v.alternate_bases)
          and h.start = v.start
          and h.reference_bases = v.reference_bases
          and h.alt_bases = TO_JSON_STRING(v.alternate_bases)
    order by v.start, v.reference_bases, TO_JSON_STRING(v.alternate_bases);
with gcount as (
  SELECT v.reference_bases, v.start, v.reference_bases, TO_JSON_STRING(v.alternate_bases) alt_bases,
           sum(if( c.genotype[OFFSET(0)] = -1,1,0)) not_called_cnt 
  FROM `genomics-public-data.1000_genomes_phase_3.variants` as v,
  UNNEST(v.call) AS c
  group by v.start, v.reference_bases, alt_bases
)
select v.not_called_cnt, count(v.not_called_cnt) as count
  from gcount as v
  group by v.not_called_cnt;

-- Query 1 - Snowflake without materialized alleles (v1)
------------------------
with filtered_snp_calls as (
    select g.chrom reference_name, g.sample_id name, concat(g.ref, '->', g.alt[0]) as mutation
        from genotypes g
        where g.chrom not in ('X', 'Y', 'MT')
              and g.ref in ('A','C','G','T')
              and g.alt[0] in ('A','C','G','T')
              and (ARRAY_SIZE(g.alt) = 1 or (ARRAY_SIZE(g.alt) = 2 and g.alt[1] = '<*>')) 
              and (strtok(g.GT, '|', 1) != '0' or strtok(g.GT, '|', 2) != '0')
              and CHARINDEX('.', g.GT) = 0
),
mutation_type_counts as (
    select reference_name, name,
           sum(iff(mutation in ('A->G', 'G->A', 'C->T', 'T->C'), 1, 0)) as transitions,
           sum(iff(mutation in ('A->C', 'C->A', 'G->T', 'T->G', 'A->T', 'T->A', 'C->G', 'G->C'), 1, 0)) as transversions
    from filtered_snp_calls
    group by reference_name, name
)
select reference_name, name, transitions, transversions, transitions/transversions as titv
    from mutation_type_counts
    where transversions > 0
    order by titv desc, name;

-- Query 2 - Snowflake without materialized alleles (v1)
------------------------
with gcount as (
  select g.pos, g.ref, g.alt, 
         sum(iff( (allele1(g.ref, g.alt, g.GT) != g.ref or allele2(g.ref, g.alt, g.GT) != g.ref) and CHARINDEX('.', g.GT) = 0,1,0)) called_non_ref,
         sum(iff(allele1(g.ref, g.alt, g.GT) != allele2(g.ref, g.alt, g.GT) and CHARINDEX('.', g.GT) = 0,1,0)) called_het,
         sum(iff(CHARINDEX('.', g.GT) = 0,1,0)) called_cnt,
         sum(iff(CHARINDEX('.', g.GT) != 0,1,0)) not_called_cnt
    from genotypes as g
    where g.chrom = '20'
    group by g.pos, g.ref, g.alt
),
homozyg_count as (
  with prep as (
    select g.pos, g.ref, g.alt, g.GT, OBJECT_CONSTRUCT(allele1(g.ref, g.alt, g.GT), count(g.GT)) homoz
        from genotypes as g
        where allele1(g.ref, g.alt, g.GT) = allele2(g.ref, g.alt, g.GT)
              and g.chrom = '20'
        group by g.pos, g.ref, g.alt, g.GT
        order by g.pos, g.ref, g.alt, g.GT
  )
 select g.pos, g.ref, g.alt, ARRAY_AGG(g.homoz) homozyg
    from prep as g
    group by g.pos, g.ref, g.alt
    order by g.pos, g.ref, g.alt
)
select v.*,
       object_construct('homozygote_count', h.homozyg,
                        'call_rate', g.called_cnt/v.info:AN*2.0,
                        'n_called', g.called_cnt,
                        'n_not_called', g.not_called_cnt,
                        'n_het', g.called_het,
                        'n_non_ref', g.called_non_ref) as variant_qc
    from variants v,
         gcount g,
         homozyg_count h
    where g.pos = v.pos
          and g.ref = v.ref
          and g.alt = v.alt
          and h.pos = v.pos
          and h.ref = v.ref
          and h.alt = v.alt
    order by v.pos, v.ref, v.alt;


-- Query 1 - Snowflake with Materialized Alleles (v2)
------------------------
with filtered_snp_calls as (
    select g.chrom reference_name, g.sample_id name, concat(g.ref, '->', g.alt[0]) as mutation
        from genotypes_V2 g
        where g.chrom not in ('X', 'Y', 'MT')
              and g.ref in ('A','C','G','T')
              and g.alt[0] in ('A','C','G','T')
              and (ARRAY_SIZE(g.alt) = 1 or (ARRAY_SIZE(g.alt) = 2 and g.alt[1] = '<*>')) 
              and (allele1 != REF or ALLELE2 != REF)
              and ALLELE1 is not NULL and ALLELE2 is not NULL
),
mutation_type_counts as (
    select reference_name, name,
           sum(iff(mutation in ('A->G', 'G->A', 'C->T', 'T->C'), 1, 0)) as transitions,
           sum(iff(mutation in ('A->C', 'C->A', 'G->T', 'T->G', 'A->T', 'T->A', 'C->G', 'G->C'), 1, 0)) as transversions
    from filtered_snp_calls
    group by reference_name, name
)
select reference_name, name, transitions, transversions, transitions/transversions as titv
    from mutation_type_counts 
    where transversions > 0 
    order by titv desc, name;

-- Query 2 - - Snowflake with Materialized Alleles (v2)
------------------------

with gcount as (
  select g.pos, g.ref, g.alt,
         sum(iff( (allele1 != g.ref or allele2 != g.ref) and ALLELE1 is not null and allele2 is not null,1,0)) called_non_ref,
         sum(iff( allele1 != allele2 and allele1 is not null and allele2 is not null,1,0)) called_het,
         sum(iff(allele1 is not null and allele2 is not null,1,0)) called_cnt,
         sum(iff(allele1 is null or allele2 is null,1,0)) not_called_cnt 
   from genotypes_V2 as g
   where g.chrom = '20'
   group by g.pos, g.ref, g.alt
),
homozyg_count as (
  with prep as (
    select g.pos, g.ref, g.alt, g.GT, OBJECT_CONSTRUCT(allele1, count(g.GT)) homoz
    from genotypes_v2 as g
    where allele1 = allele2
       and g.chrom = '20'
    group by g.pos, g.ref, g.alt, g.GT, g.allele1 
    order by g.pos, g.ref, g.alt, g.GT
)
select g.pos, g.ref, g.alt, ARRAY_AGG(g.homoz) homozyg
  from prep as g
  group by g.pos, g.ref, g.alt 
  order by g.pos, g.ref, g.alt
)
select v.*,
       object_construct('homozygote_count', h.homozyg, 
                        'call_rate', g.called_cnt/v.info:AN*2.0,
                        'n_called', g.called_cnt,
                        'n_not_called', g.not_called_cnt,
                        'n_het', g.called_het,
                        'n_non_ref', g.called_non_ref) as variant_qc
       from variants v, 
            gcount g,
            homozyg_count h 
       where g.pos = v.pos
         and g.ref = v.ref 
         and g.alt = v.alt 
         and h.pos = v.pos 
         and h.ref = v.ref 
         and h.alt = v.alt
       order by v.pos, v.ref, v.alt;
       