-- Function delivers newline-delimited data dictionary of distinct INFO
-- field data types for all files in staging for the current load
--
-- Sample usage:     Select get_info_headers_from_stg('');
--
-- Typically used as a parameter for calling PARSE_VF_INFO(...) or PARSE_VCF_INFO_FLAT(...)
--
use schema VCF;

create or replace function GET_INFO_HEADERS_FROM_STG() 
returns varchar as
$$
    select listagg(distinct text, '\n')
    from VCF_RAW_STG
    where text like '##INFO%'
$$
;