-- Function delivers newline-delimited list of header values defined 
-- for the file named in the parameter
--
-- Sample usage:     Select get_info_headers('chr1');
--
use schema GENOME1K;

create or replace function get_info_headers (vcfnm varchar) 
returns varchar as
$$
    select array_to_string(array_agg(text), '\n')
    from headers
    where vcfname = VCFNM
    and text like '##INFO%'
$$
;

-- Function delivers newline-delimited list of header values defined 
-- for al files represented in the headers table
--
-- Sample usage:     Select get_info_headers('');
--
create or replace function get_info_headers () 
returns varchar as
$$
    select array_to_string(array_agg(text), '\n')
    from ( 
        select distinct text as text 
        from headers
        where text like '##INFO%'
    ) x
$$
;