-- Workflow for deriving a multisample VCF from overlaps
--
-- Step 1:  Create a table representing the VARIANTS derived from non-NON_REF positions in any sample.
--    Include counts for all flavors of the genotype at the position
--
--     On a 4XL with 100K Exomes, this should take < 10 mins
-- 

-- Determine if ID is really necessary as it adds cost to the query

create or replace transient table MY_VARIANTS as 
	select chrom, pos, ref, alt, gt, count(*) as cnt
	from GENOTYPES_BY_SAMPLE
	where ALT_STRING != '<NON_REF>'
    --and filter = 'PASS'    -- Remove for sample data because it is too selective (an artifact of data generation)
	group by 1,2,3,4,5
    order by 1 asc, 2 asc, 6 desc;

-- Create a filter version of the table serving the purpose of a "BED" file -- positions to extract in a multi-sample VCF.
-- In the process, construct a new comprehensive ALT value per-position based on all variants seen in the data.
create or replace transient table MY_VARIANT_LOCI as 
    select chrom, pos, ref, alt, array_to_string(alt, ',') as alt_string
    from (
      select chrom, pos, ref, array_agg(distinct f.value) as alt
      from MY_VARIANTS v,
          lateral flatten(input => v.alt) f
      group by 1,2,3
      order by 1,2
    )x
;


-- detail results for VCF - LIST using GENOKEY approach
create or replace table export_detail_calls_list as
select f.chrom, f.pos, f.ref, f.alt_string, 
  rebase_gt(g.alt, f.alt, g.gt) as gt,
  count(*) as cnt,
  listagg(sample_id,',')  WITHIN GROUP (order by sample_id) as samples
from GENOTYPES_BY_SAMPLE g, MY_VARIANT_LOCI f
where  genokey(g.chrom, g.pos) <= genokey(f.chrom, f.pos)
   and genokey(g.chrom, g.end) >= genokey(f.chrom, f.pos)
group by 1,2,3,4,5
order by 1,2
;

-- OR:  Faster do do Chromosome-at-a-time 
-- avoiding expensive expression calculation for genokey
-- until multi-chrom pos range comparisons are supported SNOW-74571

create or replace procedure compute_variant_counts (RES_NAME VARCHAR) 
returns string
language javascript
as
$$

var resTable = 'spresult_' + RES_NAME

var sqlString = `
    create or replace table ` + resTable + ` (
        chrom       varchar(2),
        pos         integer,
        ref         varchar(256),
        alt_string  varchar,
        gt          varchar(5),
        cnt         integer
        )
    ;`;

var stmt = snowflake.createStatement( {sqlText: sqlString} );
stmt.execute();

var loadedChroms = '';
var chroms = ['1', '2', '3', '4', '5','6','7','8','9','10',      
              '11','12','13','14','15','16','17','18','19','20', 
              '21','22','X','Y']

for (i in chroms) {
    sqlString = `
        insert into ` + resTable + `
        with myGenotypes as (
        select  chrom, pos, end, alt, gt, sample_id
        from genotypes_by_sample
        where chrom = '` + chroms[i] +`'
        ),
        myPositions as (
        select chrom, pos, ref, alt, alt_string 
        from my_variant_loci
        where chrom = '` + chroms[i] +`'
        )
        select  b.chrom, b.pos, b.ref, b.alt_string, 
        rebase_gt(g.alt, b.alt, g.gt) as gt, 
        -- listagg(sample_id,',')  WITHIN GROUP (order by sample_id) as samples,  -- add if generating sample_ID output
        count(*) as cnt
        from myGenotypes g, myPositions b
        where g.pos <= b.pos
        and g.end >= b.pos
        group by 1,2,3,4,5
        order by 1,2;`;

    stmt = snowflake.createStatement( {sqlText: sqlString} );
    stmt.execute();
    loadedChroms += chroms[i] + ' ';
}

return loadedChroms

$$
;
