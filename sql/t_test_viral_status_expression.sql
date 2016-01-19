#Compute the t-test statistic for the differential expression of viral infected vs non viral infected 
SELECT 
p.gene as gene, 
p.study as study, 
ABS(p.x - o.y) as mean_diff, 
p.x  as x, 
p.sx2 as sx2, 
p.nx as nx, 
o.y as y, 
o.sy2 as sy2, 
o.ny as ny, 
(p.x-o.y) / SQRT((p.sx2/p.nx) + (o.sy2/o.ny)) as T 
FROM ( 
  SELECT 
  Study as study, 
  HGNC_gene_symbol as gene, 
  AVG(LOG2(normalized_count+1)) as x, 
  POW(STDDEV(LOG2(normalized_count+1)),2) as sx2, 
  COUNT(ParticipantBarcode) as nx 
  FROM 
  [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM] 
  WHERE 
  SampleTypeLetterCode = 'TP' 
  // AND Study = '_QUERY_STUDY_' 
  AND (ParticipantBarcode IN ( 
    SELECT ParticipantBarcode 
    FROM [isb-cgc:tcga_201510_alpha.Clinical_data] 
    WHERE hpv_status = 'Positive'
    GROUP BY ParticipantBarcode
  )) /* end table of participant */ 
    GROUP BY gene, study 
) as p 
JOIN ( 
  SELECT 
  Study as study, 
  HGNC_gene_symbol as gene, 
  AVG(LOG2(normalized_count+1)) AS y, 
  POW(STDDEV(LOG2(normalized_count+1)),2) AS sy2, 
  COUNT(ParticipantBarcode) as ny 
  FROM 
  [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM] 
  WHERE 
  SampleTypeLetterCode = 'TP' 
  //AND Study = '_QUERY_STUDY_' 
  AND (ParticipantBarcode IN ( 
      SELECT ParticipantBarcode 
      FROM [isb-cgc:tcga_201510_alpha.Clinical_data] 
      WHERE hpv_status = 'Negative'
      GROUP BY ParticipantBarcode
  )) /* end getting table of participants */ 
    GROUP BY 
  study, gene 
) AS o 
ON 
p.gene = o.gene and p.study = o.study 
GROUP BY gene, study, x, sx2, nx, y, sy2, ny, T, mean_diff 