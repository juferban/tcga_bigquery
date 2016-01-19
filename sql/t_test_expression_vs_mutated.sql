#Compute the t-test statistic for the differential expression of mutated vs non mutated samples 
SELECT 
   p.mut_gene as gene, 
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
  e.Study as study, 
  e.HGNC_gene_symbol as gene, 
  m.Hugo_symbol as mut_gene,
  AVG(LOG2(e.normalized_count+1)) as x, 
  POW(STDDEV(LOG2(e.normalized_count+1)),2) as sx2, 
  COUNT(e.ParticipantBarcode) as nx 
FROM (
/* New code to get expression for all mutated genes and samples */
  SELECT 
    HGNC_gene_symbol,
    Study,
    normalized_count,
    ParticipantBarcode
  FROM [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM] 
  WHERE 
    SampleTypeLetterCode = 'TP'
    AND HGNC_gene_symbol = '_QUERY_GENE_'
    //AND Study = 'BRCA' 
) e
INNER JOIN (
  SELECT 
    ParticipantBarcode, 
    Hugo_symbol
  FROM [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] as m 
  WHERE 
    Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
    //AND m.Study = 'BRCA'
    GROUP BY 
      ParticipantBarcode, Hugo_symbol 
) m /* end table of participant */ 
ON 
e.ParticipantBarcode = m.ParticipantBarcode
GROUP BY gene, mut_gene, study 
/* end of new code */
     ) as p 
     
 JOIN ( 
    SELECT 
  ee.Study as study, 
  ee.gene as gene, 
  ee.mut_gene as mut_gene,
  AVG(LOG2(ee.normalized_count+1)) as y, 
  POW(STDDEV(LOG2(ee.normalized_count+1)),2) as sy2, 
  COUNT(ee.ParticipantBarcode) as ny
FROM (
// Cartesian join to get all gene pair combinations for each study
  SELECT 
    em.Study as Study, 
    em.gene as gene, 
    gm.Hugo_symbol as mut_gene,
    em.ParticipantBarcode as ParticipantBarcode,
    em.normalized_count as normalized_count
  FROM (  
// First we need to filter the expression table to get only those samples that have been analyzed for mutation
    SELECT 
      e.Study as Study, 
      e.HGNC_gene_symbol as gene, 
      e.ParticipantBarcode as ParticipantBarcode,
      e.normalized_count as normalized_count
    FROM (
      SELECT 
        HGNC_gene_symbol,
        Study,
        normalized_count,
        ParticipantBarcode
      FROM [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM] 
      WHERE 
        SampleTypeLetterCode = 'TP'
        AND HGNC_gene_symbol = '_QUERY_GENE_'
        //AND Study = 'BRCA' 
    ) e
    JOIN (
      SELECT ParticipantBarcode
      FROM [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls]
      GROUP BY ParticipantBarcode
    ) ms
    ON e.ParticipantBarcode = ms.ParticipantBarcode
  ) em
  JOIN (
#Get list of mutated genes per study
    SELECT Hugo_symbol, Study
    FROM [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls]
    WHERE Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
      //AND Study = 'BRCA'
    GROUP by Hugo_symbol, Study
  ) gm
  ON em.Study = gm.Study
) ee
LEFT OUTER JOIN (
// Get samples that have mutations for each gene and do an "Anti join to remove those samples
  SELECT 
    ParticipantBarcode, 
    Hugo_symbol
  FROM [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] as m 
  WHERE 
    Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
    //AND m.Study = 'BRCA'
    GROUP BY 
      ParticipantBarcode, Hugo_symbol 
) m /* end table of participant */ 
ON 
ee.ParticipantBarcode = m.ParticipantBarcode
AND ee.mut_gene = m.Hugo_symbol
WHERE m.ParticipantBarcode is NULL
GROUP BY gene, mut_gene, study 
  ) AS o 
 ON 
   p.mut_gene = o.mut_gene and p.study = o.study 
 GROUP BY gene, study, x, sx2, nx, y, sy2, ny, T, mean_diff 