-- Scalar UDF to convert the INFO field from VCF entries WITH SINGLETON ALT VALUES PER POSITION
-- ****************************************************************************************************
-- THIS DIFFERES FROM THE VERSION USED IN 1000Genomes SCHEMA, WHICH HANDLES MULTIPLE ALT VALUES PER POS
-- ****************************************************************************************************
-- 
-- Converts INFO into a JSON representation for easier querying (or to prepare for column extraction)
-- 
-- This version is appropriate for some annotation VCFs that never have multiple-alt alleles such as GNOMAD
--
-- Inputs:
--     info ==>  the INFO value or column of a VCF table, which is semicolon-separated, fields defined by '=' -- in format such as 'RS=1052373574;dbSNPBuildID=150;SSR=0;...'
--     struct ==>  rows of text (newline separated) containing the ##INFO header entries (from the original VCF file) containing metadata for the INFO field
--         Use the helper function get_info_header to extract this text -- e.g.   get_info_header('my_vcf_name')
--  So example calling would be 
--      select parse_vcf_info_flat('AC=2130;AF=0.425319;AN=5008;NS=2504;DP=103152;EAS_AF=0.3363;AMR_AF=0.3602;AFR_AF=0.4909;EUR_AF=0.4056;SAS_AF=0.4949;AA=|||unknown(NO_COVERAGE);VT=INDEL'
--                          , get_info_headers('chr1'));
--
--      or applied to an entire table:
--      select parse_vcf_info_flat(info, get_info_headers_from_stg()) from VCF_RAW_STG ;
-- *************************************************************************************************************

create or replace function PARSE_VCF_INFO_FLAT(info varchar, struct varchar)
  RETURNS variant
  LANGUAGE JAVASCRIPT
  AS '
    var setup = function() {
        structMap = new Map();
        infos = [];
        
        parsedValue = (theKey, theValue) => {
            return  structMap[theKey].type == "Integer" ? Number.parseInt(theValue) :
                    structMap[theKey].type == "Float"   ? Number.parseFloat(theValue) :
                    theValue;
        };

        STRUCT.split("\\n").forEach( (row, i)=> {
            info = (row.substring(row.indexOf("<")+1, row.indexOf(">"))).split(",")
            infos[i] = {  ID : info[0].split("ID=")[1], 
              Number : info[1].split("Number=")[1], 
              Type : info[2].split("Type=")[1], 
              Description : info[3].split("Description=")[1].replace(/"/g,"")
            };
        });

        infos.forEach((elem) => {
        structMap[elem.ID] = {number: elem.Number, type: elem.Type};
        });
    }
        
    if (typeof(setup_done) === "undefined") {
        setup();
        setup_done = true;  // setting global variable to true
    }
    
    infoArray = INFO.split(";");
    infoMap = new Map();

    infoArray.forEach((item) => {
        theKey = item.split("=",2)[0].trim();

        if (structMap[theKey].number == "0") {
        theValue = true;
        }
        else if (structMap[theKey].number == "1" || structMap[theKey].number == "A") {
        theValue = parsedValue(theKey,item.split("=",2)[1].trim());
        }
        else {
        theValue = item.split("=",2)[1].trim().split(",").map(val => parsedValue(theKey, val));
        }

        infoMap[theKey] = theValue;
    });

    return(infoMap);
'
;