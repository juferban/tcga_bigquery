#'  
#'  
#' # Defining cohorts and performing differential expression analysis 
#'  
#' Differential expression (DE) is a common analysis that determines if the mean 
#' gene expression level differ between two groups. Most often we have two 
#' groups, normal and disease. But, in the case of TCGA, normal tissue samples 
#' are not always available. Nonetheless, samples can be grouped in many different 
#' ways, and with the number of samples and the types of data, one can be quite creative with 
#' analysis. In this work, I will be creating groups of samples based on somatic mutations, 
#' changes in the DNA of tumors, and associating it with differentially expressed genes. 
#'  
 
## ----set project,-------------------------------------------------------- 
#project <- YOUR PROJECT ID 

#'
#' The first question: "what types of mutation data are available?" 
#'   
library(dplyr) 
library(bigrquery) 
library(scales) 
library(ggplot2) 
library(tidyr)
#'library(ISBCGCExamples) 
 

q <- "SELECT 
         Variant_Classification, COUNT(*) AS n 
       FROM 
         [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] 
       WHERE 
         Study = 'BRCA' 
       GROUP BY 
         Variant_Classification" 
 

result <- query_exec(q, project) 
result 
 

#'  
#'  
#'  
#' This query is selecting all the variant classification labels for samples 
#' in the BRCA study using the Somatic Mutation table. The GROUP BY portion of the query 
#' ensures we only get each category once, similar to the DISTINCT keyword. 
#'  
#' Oftentimes, people are interested in specific types of variants, such as exonic variants. 
#' Here, I'm interested in whether any protein altering variant might be associated with 
#' differential expression. In order to perform such an analysis, we 
#' need to define two groups. For now, let's focus on the PIK3CA gene. So
#' let's find out how many individuals have a mutation in PIK3CA, in the BRCA study. 
#'  
## ----Get list of samples with mutations in PIK3CA ----------------------------------------------------- 
q <- "SELECT 
         ParticipantBarcode 
       FROM 
         [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] 
       WHERE 
         Hugo_Symbol = 'PIK3CA' 
         AND Study = 'BRCA' 
         AND Variant_Classification not in (\"Intron\",\"RNA\",\"IGR\",\"lincRNA\",\"3'UTR\",\"Silent\",\"5'UTR\")
       GROUP BY 
         ParticipantBarcode" 
 

query_exec(q, project) 

# The directory in which the files containing SQL reside.
sqlDir = file.path("Path/to/sql/files/")
 
#'  
#' In this query, we get that a fair number of Samples in BRCA are mutated for PIK3CA which make this gene
#' a suitable candidate for our test
#' 
#' Now we are going to compare the expression for all genes in the mutated samples vs the samples that don't
#' have a mutation for the query gene. We are going to group the data by gene and study (cancer type)
#'  
## ----chunk5-------------------------------------------------------------- 
ptm1 <- proc.time() 

query_gene <- "PIK3CA"

# Now we are ready to run the t-test query query. 
result = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "t_test_mut_vs_expression.sql"), 
  project=project, 
  replacements=list("_QUERY_GENE_"= query_gene ) 
                    
) 

ptm2 <- proc.time() - ptm1 
cat("Wall-clock time for BigQuery:",ptm2[3]) 


compute_df <- function(d) { 
 ((d$sx2/d$nx + d$sy2/d$ny)^2) / 
 ((1/(d$nx-1))*(d$sx2/d$nx)^2 + (1/(d$ny-1))*(d$sy2/d$ny)^2) 
} 

result$df <- compute_df(result) 
result$p_value <- sapply(1:nrow(result), function(i) 2*pt(abs(result$T[i]), result$df[i],lower=FALSE)) 
result$fdr <- p.adjust(result$p_value, "fdr") 
result$gene_label <- 1 
result$gene_label[result$gene == query_gene] <- 2 

result_matrix <- result %>% dplyr::select(gene, study, p_value) %>% tidyr::spread(study, p_value) 

# ordered by difference in means 
head(result) 
 

# ordered by T statistic 
head(result[order(result$T, decreasing=T), ]) 
 
# ordered by FDR 
head(result[order(result$fdr, decreasing=F), ]) 
 

qplot(data=result, x=T, y=mean_diff, shape=as.factor(gene_label), col=as.factor(fdr < 0.05), ylab="Difference in Means", xlab="T statistic") 


# Repeat the Analysis with "Optimized" code
ptm1_alt <- proc.time() 

#query_gene <- "PIK3CA"

# Now we are ready to run the t-test query query. 
result_alt = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "t_test_mut_vs_expression_alt.sql"), 
  project=project, 
  replacements=list("_QUERY_GENE_"= query_gene ) 
  
) 

ptm2_alt <- proc.time() - ptm1_alt 
cat("Wall-clock time for BigQuery:",ptm2_alt[3]) 


result_alt$df <- compute_df(result_alt) 
result_alt$p_value <- sapply(1:nrow(result_alt), function(i) 2*pt(abs(result_alt$T[i]), result_alt$df[i],lower=FALSE)) 
result_alt$fdr <- p.adjust(result_alt$p_value, "fdr") 
result_alt$gene_label <- 1 
result_alt$gene_label[result_alt$gene == query_gene] <- 2 

result_matrix_alt <- result_alt %>% dplyr::select(gene, study, p_value) %>% tidyr::spread(study, p_value) 

# ordered by difference in means 
head(result_alt) 


# ordered by T statistic 
head(result_alt[order(result_alt$T, decreasing=T), ]) 

# ordered by FDR 
head(result_alt[order(result_alt$fdr, decreasing=F), ]) 


qplot(data=result_alt, x=T, y=mean_diff, shape=as.factor(gene_label), col=as.factor(fdr < 0.05), ylab="Difference in Means", xlab="T statistic") 




#'  
#' Let's check on a single gene. We'll pull down 
#' the actual expression values, and use the R T-test. 
#'  
#'  
## ----ttest_fig2---------------------------------------------------------- 
 

q <- "SELECT HGNC_gene_symbol, ParticipantBarcode, LOG2(normalized_count+1) 
 FROM [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM] 
 WHERE Study = 'BRCA' 
 and HGNC_gene_symbol = 'PIK3CA' 
 and SampleTypeLetterCode = 'TP' 
 and (ParticipantBarcode IN ( 
   SELECT 
     m.ParticipantBarcode, 
   FROM 
     [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] as m 
   WHERE 
     m.Hugo_Symbol = 'PIK3CA' 
     AND m.Study = 'BRCA' 
   GROUP BY 
     m.ParticipantBarcode 
   )) 
 " 
 mutExpr <- query_exec(q, project)   # SOME DUPLCIATES 
 

q <- " 
 SELECT HGNC_gene_symbol, ParticipantBarcode, LOG2(normalized_count+1) 
 FROM [isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM] 
 WHERE Study = 'BRCA' 
 and HGNC_gene_symbol = 'PIK3CA' 
 and SampleTypeLetterCode = 'TP' 
 and (ParticipantBarcode IN ( 
   SELECT 
     ParticipantBarcode 
   FROM 
     [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] 
   WHERE Study = 'BRCA' 
   and (ParticipantBarcode NOT IN ( 
       SELECT 
         m.ParticipantBarcode, 
       FROM 
         [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls] as m 
       WHERE 
         m.Hugo_Symbol = 'PIK3CA' 
         and m.Study = 'BRCA' 
       GROUP BY 
         m.ParticipantBarcode 
       )) /* end getting table of participants without mutations */ 
   GROUP BY ParticipantBarcode 
 )) /* end getting table of participants */ 
 " 
 wtExpr <- query_exec(q, project)   # SOME DUPLCIATES 
 

 t.test(mutExpr$f0_, wtExpr$f0_) 
 

boxplot(list(Mutation_In_PIK3CA=mutExpr$f0_, No_Mutation_In_PIK3CA=wtExpr$f0_), ylab="PIK3CA LOG2 expression") 

#' Know let's investigate the relation between expression in one gene vs the mutation status
#' of all genes to try to identify new mutations that may affect the expression of the query gene
#' 

ptm3 <- proc.time() 

query_gene <- "PIK3CA"

# Now we are ready to run the t-test query query. 
result0 = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "t_test_expression_vs_mutated.sql"), 
  project=project, 
  replacements=list("_QUERY_GENE_"= query_gene ) 
  
) 

ptm4 <- proc.time() - ptm3 
cat("Wall-clock time for BigQuery:",ptm4[3]) 


compute_df <- function(d) { 
  ((d$sx2/d$nx + d$sy2/d$ny)^2) / 
    ((1/(d$nx-1))*(d$sx2/d$nx)^2 + (1/(d$ny-1))*(d$sy2/d$ny)^2) 
} 

result0$df <- compute_df(result0) 
result0$p_value <- sapply(1:nrow(result0), function(i) 2*pt(abs(result0$T[i]), result0$df[i],lower=FALSE)) 
result0$fdr <- p.adjust(result0$p_value, "fdr") 
result0$gene_label <- 1 
result0$gene_label[result0$gene == query_gene] <- 2 

result_matrix0 <- result0 %>% dplyr::select(gene, study, p_value) %>% tidyr::spread(study, p_value) 

# ordered by difference in means 
head(result) 


# ordered by T statistic 
head(result[order(result$T, decreasing=T), ]) 

# ordered by FDR 
head(result[order(result$fdr, decreasing=F), ]) 


qplot(data=result, x=T, y=mean_diff, shape=as.factor(gene_label), col=as.factor(fdr < 0.05), ylab="Difference in Means", xlab="T statistic") 

# Repeat with "Optimized" code

ptm3_alt <- proc.time() 

#query_gene <- "PIK3CA"

# Now we are ready to run the t-test query query. 
result0_alt = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "t_test_expression_vs_mutated_alt.sql"), 
  project=project, 
  replacements=list("_QUERY_GENE_"= query_gene ) 
  
) 

ptm4_alt <- proc.time() - ptm3_alt 
cat("Wall-clock time for BigQuery:",ptm4_alt[3]) 


result0_alt$df <- compute_df(result0_alt) 
result0_alt$p_value <- sapply(1:nrow(result0_alt), function(i) 2*pt(abs(result0_alt$T[i]), result0_alt$df[i],lower=FALSE)) 
result0_alt$fdr <- p.adjust(result0_alt$p_value, "fdr") 
result0_alt$gene_label <- 1 
result0_alt$gene_label[result0_alt$gene == query_gene] <- 2 

result_matrix0_alt <- result0_alt %>% dplyr::select(gene, study, p_value) %>% tidyr::spread(study, p_value) 

# ordered by difference in means 
head(result0_alt) 


# ordered by T statistic 
head(result0_alt[order(result0_alt$T, decreasing=T), ]) 

# ordered by FDR 
head(result0_alt[order(result0_alt$fdr, decreasing=F), ]) 


qplot(data=result0_alt, x=T, y=mean_diff, shape=as.factor(gene_label), col=as.factor(fdr < 0.05), ylab="Difference in Means", xlab="T statistic") 



#'  
#' Let's get count for gene pair mutation status  
#' and perform a fisher test to evaluate comutated or mutually exclusive mutated genes. 
#'  
#'  
## ----chunk5-------------------------------------------------------------- 
ptm5 <- proc.time() 

query_gene <- "PIK3CA"

# Now we are ready to run the t-test query query. 
result1 = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "mut_vs_mut_counts.sql"), 
  project=project, 
  replacements=list("_QUERY_GENE_"= query_gene ) 
) 

ptm6 <- proc.time() - ptm5 
cat("Wall-clock time for BigQuery:",ptm6[3]) 

#result1 <- result1 %>% dplyr::filter(gene != "Unknown" & mm_count > 5 & no_mm_count >5 & m_nom_count >5 & nono_count >5)

system.time(result1$p_value <- sapply(1:nrow(result1), function(i) fisher.test(matrix(c(result1$mm_count[i], 
                                                                            result1$no_mm_count[i],
                                                                            result1$m_nom_count[i], 
                                                                            result1$nono_count[i]), 
                                                                            nrow=2), 
                                                                          alternative = "greater")$p.value))



result_matrix1 <- result1 %>% dplyr::select(gene, study, p_value) %>% tidyr::spread(study, p_value) 


sessionInfo() 


