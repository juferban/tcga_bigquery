# Compute the correlation between expression data. 
select query_gene, 
       target_gene, 
       SampleTypeLetterCode, 
       cancer_type, 
       count(distinct(SampleBarcode)) AS num_observations, 
       corr(log2_count2, log2_count) AS pearson_corr
from (
  SELECT SampleBarcode, 
         SampleTypeLetterCode, 
         AliquotBarcode,
         Study                as cancer_type, 
         HGNC_gene_symbol     as target_gene, 
         log2(normalized_count+1) As log2_count,
         case 
           when HGNC_gene_symbol = '_QUERY_GENE_' then 0 
           else 1 
         end as flag,
         '_QUERY_GENE_' as query_gene,
         first_value(log2(normalized_count + 1)) over (partition by SampleBarcode, 
                                              SampleTypeLetterCode,        
                                              AliquotBarcode 
                                     order by flag, log2_count) log2_count2
  FROM   [_EXPRESSION_TABLE_] 
) 
where target_gene is not NULL
group by query_gene, target_gene, cancer_type, SampleTypeLetterCode
# HAVING 
#   # although we will compute correlations for all genes pairs with EGFR, we only want to keep values where 
#   # the correlation was based on at least 30 observations 
#   num_observations >= _MINIMUM_NUMBER_OF_OBSERVATIONS_ 
order by query_gene, target_gene, cancer_type, SampleTypeLetterCode

