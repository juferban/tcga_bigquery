#Compute the fisher-test statistic for comparing mutation status of gene pairs. 
SELECT 
  t.study as study,
  c.hugo_symbol as gene,
  c.mm_count as mm_count,
  c.no_mm_count as no_mm_count,
  (m.mut_samples - mm_count) as m_nom_count,
  (t.sample_count - (mm_count + no_mm_count + (m.mut_samples - mm_count))) as nono_count,
  t.sample_count as sample_count
FROM (
  SELECT 
    count(distinct ParticipantBarcode) as sample_count, 
    Study
  FROM 
    [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] 
  GROUP BY 
    Study 
) t
JOIN(
  SELECT
    mm.study as study,
    mm.hugo_symbol as hugo_symbol,
    mm.mm_count as mm_count,
    no_mm.no_mm_count as no_mm_count
  FROM (
    SELECT 
      Study as study, 
      HUGO_symbol,  
      COUNT(distinct ParticipantBarcode) as mm_count
    FROM 
      [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] 
    WHERE 
      Tumor_SampleTypeLetterCode = 'TP' 
      //AND Study = 'BRCA' 
      AND Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
      AND (ParticipantBarcode IN ( 
        SELECT 
          m.ParticipantBarcode, 
        FROM 
          [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] as m 
        WHERE 
          m.Hugo_Symbol = '_QUERY_GENE_'
          AND Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
          //AND m.Study = 'BRCA' 
        GROUP BY 
          m.ParticipantBarcode 
        )
      ) /* end table of participant */ 
    GROUP BY hugo_symbol, study 
  ) as mm
  JOIN (
  /* get samples count for genes that have mutation in the target gene but not in the query gene */
    SELECT 
      Study as study, 
      HUGO_symbol,  
      COUNT(distinct ParticipantBarcode) as no_mm_count 
    FROM 
      [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] 
    WHERE 
      Tumor_SampleTypeLetterCode = 'TP' 
      //AND Study = 'BRCA' 
      AND Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
      AND (ParticipantBarcode not IN ( 
        SELECT 
          m.ParticipantBarcode, 
        FROM 
          [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] as m 
        WHERE 
          m.Hugo_Symbol = '_QUERY_GENE_'
          AND Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
          //AND m.Study = 'BRCA' 
        GROUP BY 
          m.ParticipantBarcode 
        )
      ) /* end table of participant */ 
    GROUP BY hugo_symbol, study 
  ) no_mm
  ON 
  mm.hugo_symbol = no_mm.hugo_symbol and mm.study = no_mm.study 
) c
ON t.study = c.study
JOIN (
  SELECT 
    count(distinct ParticipantBarcode) as mut_samples, Study
  FROM 
    [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] as m 
  WHERE 
    m.Hugo_Symbol = '_QUERY_GENE_'
    AND Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
    GROUP BY 
    Study 
) m
ON
t.study = m.study