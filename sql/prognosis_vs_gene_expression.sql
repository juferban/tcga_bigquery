/* Compute the correlation between expression and prognosis data (measure as survival days from diagnosis). */
select a_exp.target_gene, 
       a_exp.SampleTypeLetterCode, 
       a_exp.cancer_type, 
       count(distinct(a_exp.SampleBarcode)) AS num_observations, 
       corr(a_exp.log2_count, b_clin.days_to_last_known_alive) AS pearson_corr
from (
  SELECT ParticipantBarcode,
         SampleBarcode, 
         SampleTypeLetterCode, 
         AliquotBarcode,
         Study                as cancer_type, 
         HGNC_gene_symbol     as target_gene, 
         log2(normalized_count+1) As log2_count,         
  FROM   [_EXPRESSION_TABLE_]
) as a_exp
join (
  SELECT ParticipantBarcode,
         days_to_last_known_alive 
  FROM [_CLINICAL_TABLE_] 
) as b_clin
on a_exp.ParticipantBarcode = b_clin.ParticipantBarcode
where a_exp.target_gene is not NULL
group by a_exp.target_gene, a_exp.cancer_type, a_exp.SampleTypeLetterCode
