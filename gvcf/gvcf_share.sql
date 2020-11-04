use role ACCOUNTADMIN;
use schema GVCF;

---To create share "GVCF" with relevant tables and views

CREATE SHARE "GVCF" COMMENT='';
GRANT USAGE ON DATABASE "GENOME" TO SHARE "GVCF";
GRANT USAGE ON SCHEMA "GENOME"."GVCF" TO SHARE "GVCF";
GRANT SELECT ON VIEW "GENOME"."GVCF"."CLINVAR" TO SHARE "GVCF";
GRANT SELECT ON VIEW "GENOME"."GVCF"."GENOTYPES_BY_POSITION" TO SHARE "GVCF";
GRANT SELECT ON VIEW "GENOME"."GVCF"."GENOTYPES_BY_SAMPLE" TO SHARE "GVCF";
GRANT SELECT ON VIEW "GENOME"."GVCF"."GENOTYPES_ALT" TO SHARE "GVCF";
GRANT SELECT ON VIEW "GENOME"."GVCF"."HEADERS" TO SHARE "GVCF";
GRANT SELECT ON VIEW "GENOME"."GVCF"."SNP" TO SHARE "GVCF";



----To add UDFs to "GVCF" share
GRANT USAGE ON FUNCTION allele1(varchar, variant, varchar)
to share "GVCF";
GRANT USAGE ON FUNCTION allele2(varchar, variant, varchar)
to share "GVCF";
GRANT USAGE ON FUNCTION  allele_position(varchar, int)
to share "GVCF";
GRANT USAGE ON FUNCTION  is_homozygous(varchar)
to share "GVCF";
GRANT USAGE ON FUNCTION  at_location(integer, integer, integer)
to share "GVCF";


----To share with specific account name
ALTER SHARE "GVCF" ADD ACCOUNTS = ***AccountName***;
