/* Alternative implementation with similar performance */ 
  select study, 
mut_gene, 
ABS(x - y) as mean_diff,
x, 
sx2, 
nx, 
y, 
sy2, 
ny, 
(x-y) / SQRT((sx2/nx) + (sy2/ny)) as T 
from (
  -- pivot rows to columns
  select study, mut_gene, 
  max(if (mtype = 'mutated'    , a,   null)) as x   ,
  max(if (mtype = 'mutated'    , sa2, null)) as sx2 ,
  max(if (mtype = 'mutated'    , na,  null)) as nx  ,
  max(if (mtype = 'not mutated', a,   null)) as y   ,
  max(if (mtype = 'not mutated', sa2, null)) as sy2 ,
  max(if (mtype = 'not mutated', na,  null)) as ny
  from   (
    SELECT  ee.Study                                   as study, 
    ee.gene                                    as gene, 
    ee.mut_gene                                as mut_gene,
    if(m.hugo_symbol is null, 'not mutated','mutated') mtype,
    AVG(LOG2(ee.normalized_count+1))           as a, 
    POW(STDDEV(LOG2(ee.normalized_count+1)),2) as sa2, 
    COUNT(ee.ParticipantBarcode)               as na
    FROM   (
      -- Cartesian join to get all gene pair combinations for each study
      SELECT em.Study              as Study, 
      em.gene               as gene, 
      gm.Hugo_symbol        as mut_gene,
      em.ParticipantBarcode as ParticipantBarcode,
      em.normalized_count   as normalized_count
      FROM   (  
        -- First we need to filter the expression table to get only those samples that have been analyzed for mutation
        SELECT e.Study              as Study, 
        e.HGNC_gene_symbol   as gene, 
        e.ParticipantBarcode as ParticipantBarcode,
        e.normalized_count   as normalized_count
        FROM (
          SELECT HGNC_gene_symbol,
          Study,
          normalized_count,
          ParticipantBarcode
          FROM   [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM] 
          WHERE  SampleTypeLetterCode = 'TP'
          AND    HGNC_gene_symbol = '_QUERY_GENE_'
        ) e
        JOIN (
          SELECT ParticipantBarcode
          FROM   [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls]
          GROUP  BY ParticipantBarcode
        ) ms
        ON e.ParticipantBarcode = ms.ParticipantBarcode
      ) em
      JOIN   (
        -- Get list of mutated genes per study
        SELECT Hugo_symbol, Study
        FROM [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls]
        WHERE Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
        GROUP by Hugo_symbol, Study
      ) gm
      ON em.Study = gm.Study                                                                      
    ) ee
    LEFT OUTER JOIN 
    -- join actual mutations so that you can determine actual mutations versus non mutated genes
    (                 
    SELECT ParticipantBarcode, 
    Hugo_symbol
    FROM   [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] as m 
    WHERE  Variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
    GROUP  BY ParticipantBarcode, Hugo_symbol 
    ) m /* end table of participant */ 
      ON    ee.ParticipantBarcode = m.ParticipantBarcode AND ee.mut_gene = m.Hugo_symbol
    GROUP BY gene, mut_gene, study, mtype
  )      
  group by study, mut_gene       
)
where x is not null and y is not null
order by 1,2
