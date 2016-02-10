/* Alternative implementation  with slightly higher performance */ 
/* This version is not filtering out samples that have expression that have NOT been analyzed for mutations. */
select study, hgnc_gene_symbol as gene,
abs(x-y) as mean_diff,
x,
sx2,
nx,
y,
sy2,
ny,
(x-y) / SQRT((sx2/nx) + (sy2/ny)) as T
from   (
  select study, hgnc_gene_symbol,
  max(case when query_gene = 1 then avg else null end) x   ,
  max(case when query_gene = 1 then std else null end) sx2 ,
  max(case when query_gene = 1 then cnt else null end) nx  ,
  max(case when query_gene = 0 then avg else null end) y   ,
  max(case when query_gene = 0 then std else null end) sy2 ,
  max(case when query_gene = 0 then cnt else null end) ny
  from   (select r.study,
          r.HGNC_gene_symbol,
          query_gene,
          AVG(LOG2(normalized_count+1))           as avg,
          POW(STDDEV(LOG2(normalized_count+1)),2) as std,
          COUNT(r.ParticipantBarcode)             as cnt       
          from   [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM]     r
          join   (
            SELECT ParticipantBarcode,
            case when hugo_symbol = '_QUERY_GENE_' then 1 else 0 end query_gene,
            FROM   [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls]        
            WHERE   variant_classification not in ("Intron","RNA","IGR","lincRNA","3'UTR","Silent","5'UTR")
            AND SampleTypeLetterCode = 'TP'
            GROUP  BY ParticipantBarcode, query_gene
          ) m
          on     r.ParticipantBarcode = m.ParticipantBarcode
          where SampleTypeLetterCode = 'TP'
          group by r.study,
          r.HGNC_gene_symbol,
          query_gene
  )
  group by study, hgnc_gene_symbol  having nx is not null and ny is not null    
)        
;

