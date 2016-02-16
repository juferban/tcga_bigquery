SELECT 
   p.target_gene as gene, 
   p.cancer_type as cancer_type, 
   ABS(p.x - t.y) as mean_diff, 
   p.x  as x, 
   p.sx2 as sx2, 
   p.nx as nx, 
   t.y as y, 
   t.sy2 as sy2, 
   t.ny as ny, 
   (p.x-t.y) / SQRT((p.sx2/p.nx) + (t.sy2/t.ny)) as T 
FROM ( 
  SELECT a_exp.target_gene, 
       a_exp.SampleTypeLetterCode, 
       a_exp.cancer_type, 
       avg(log2_count) as x,
       pow(stddev(log2_count),2) as sx2,
       count(distinct(a_exp.SampleBarcode)) AS nx, 
  FROM (
    SELECT ParticipantBarcode,
           SampleBarcode, 
           SampleTypeLetterCode, 
           AliquotBarcode,
           Study                as cancer_type, 
           HGNC_gene_symbol     as target_gene, 
           log2(normalized_count+1) As log2_count,         
    FROM   [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
  ) as a_exp
  JOIN (
    SELECT ParticipantBarcode,
           days_to_last_known_alive,
    FROM [isb-cgc:tcga_201510_alpha.Clinical_data] 
    WHERE gender = 'MALE'
  ) as b_clin
  on a_exp.ParticipantBarcode = b_clin.ParticipantBarcode
  where a_exp.target_gene is NOT NULL
  group by a_exp.target_gene, a_exp.cancer_type, a_exp.SampleTypeLetterCode
) as p     
JOIN (
 SELECT a_exp.target_gene, 
       a_exp.SampleTypeLetterCode, 
       a_exp.cancer_type, 
       avg(log2_count) as y,
       pow(stddev(log2_count),2) as sy2,
       count(distinct(a_exp.SampleBarcode)) AS ny, 
  FROM (
    SELECT ParticipantBarcode,
         SampleBarcode, 
         SampleTypeLetterCode, 
         AliquotBarcode,
         Study                as cancer_type, 
         HGNC_gene_symbol     as target_gene, 
         log2(normalized_count+1) As log2_count,         
    FROM   [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]
  ) as a_exp
  join (
    SELECT ParticipantBarcode,
         days_to_last_known_alive,
    FROM [isb-cgc:tcga_201510_alpha.Clinical_data] 
    WHERE gender = 'FEMALE'
  ) as b_clin
  on a_exp.ParticipantBarcode = b_clin.ParticipantBarcode
  where a_exp.target_gene is not NULL
  group by a_exp.target_gene, a_exp.cancer_type, a_exp.SampleTypeLetterCode
) AS t 
ON 
p.target_gene = t.target_gene and p.cancer_type = t.cancer_type 
GROUP BY gene, cancer_type, x, sx2, nx, y, sy2, ny, T, mean_diff