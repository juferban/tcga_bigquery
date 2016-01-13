# This query performs a Spearman (ranked) correlation between the RNA-Seq expression in one gene - and the RNA-expression data 
# partitioned based on gene pair, sample type and cancer type.  The correlation is done by the BigQuery CORR function, but first 
# we need to use the RANK function to turn expression values into ranks. 

#### QUESTION: should the COUNT(1) be COUNT(2) instead ??? 

 
SELECT 
   query_gene, 
   target_gene, 
   SampleTypeLetterCode, 
   cancer_type, 
   count(distinct(SampleBarcode)) AS num_observations, 
   CORR(query_expr_rank, target_expr_rank) AS spearman_corr 
FROM ( 
# in order to do a rank-based (Spearman) correlation, we need to turn the expression values into ranks 
  SELECT 
     SampleBarcode,
     SampleTypeLetterCode,
     AliquotBarcode,
     cancer_type,
     query_gene, 
     target_gene, 
     RANK() OVER (PARTITION BY SampleTypeLetterCode, cancer_type, query_gene, target_gene ORDER BY log2_count2 ASC) AS query_expr_rank, 
     RANK() OVER (PARTITION BY SampleTypeLetterCode, cancer_type, query_gene, target_gene ORDER BY log2_count ASC) AS target_expr_rank, 
  FROM (
    SELECT   
      SampleBarcode, 
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
)
GROUP BY 
  query_gene, 
  target_gene,
  SampleTypeLetterCode, 
  cancer_type
# HAVING 
#   # although we will compute correlations for all genes and proteins, we only want to keep values where 
#   # the correlation was based on at least 30 observations 
#   num_observations >= _MINIMUM_NUMBER_OF_OBSERVATIONS_ 
ORDER BY 
   # and finally we want to order the ouputs from largest positive correlation, descending to the most negative 
   # correlation 
spearman_corr DESC 
