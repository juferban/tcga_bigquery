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

# The directory in which the files containing SQL reside.
sqlDir = file.path("/Path/to/sql/files")

#'  
#' In this query, 
#' 
#' Now we are going to compare the expression for all genes in the mutated samples vs the samples that don't
#' have a mutation for the query gene. We are going to group the data by gene and study (cancer type)
#'  
## ----chunk5-------------------------------------------------------------- 
expressionTable <- 'isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM'
clinicalTable <- 'isb-cgc:tcga_201510_alpha.Clinical_data'

ptm1 <- proc.time() 


# Now we are ready to run the t-test query query. 
result = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "prognosis_vs_gene_expression.sql"), 
  project=project,
  replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                    "_CLINICAL_TABLE_"=clinicalTable
  )

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

# Analysis of differentially expressed genes based on gender
ptm3 <- proc.time() 


# Now we are ready to run the t-test query query. 
result_g = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "prognosis_vs_gene_expression_gender.sql"), 
  project=project,
  replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                    "_CLINICAL_TABLE_"=clinicalTable
  )
  
) 

ptm4 <- proc.time() - ptm3 
cat("Wall-clock time for BigQuery:",ptm4[3]) 

result_g$df <- compute_df(result_g) 
result_g$p_value <- sapply(1:nrow(result_g), function(i) 2*pt(abs(result_g$T[i]), result_g$df[i],lower=FALSE)) 
result_g$fdr <- p.adjust(result_g$p_value, "fdr") 
result_g$gene_label <- 1 
result_g$gene_label[result_g$gene == query_gene] <- 2 

result_g_matrix <- result_g %>% dplyr::select(gene, cancer_type, p_value) %>% tidyr::spread(cancer_type, p_value) 

# ordered by difference in means 
head(result_g) 


# ordered by T statistic 
head(result_g[order(result_g$T, decreasing=T), ]) 

# ordered by FDR 
head(result_g[order(result_g$fdr, decreasing=F), ]) 


qplot(data=result_g, x=T, y=mean_diff, shape=as.factor(gene_label), col=as.factor(fdr < 0.05), ylab="Difference in Means", xlab="T statistic") 


# Analysis of differentially expressed genes based on treatment_status
ptm5 <- proc.time() 


# Now we are ready to run the t-test query query. 
result_t = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "prognosis_vs_gene_expression_treated.sql"), 
  project=project,
  replacements=list("_EXPRESSION_TABLE_"=expressionTable,
                    "_CLINICAL_TABLE_"=clinicalTable
  )
  
) 

ptm6 <- proc.time() - ptm5 
cat("Wall-clock time for BigQuery:",ptm6[3]) 

result_t$df <- compute_df(result_t) 
result_t$p_value <- sapply(1:nrow(result_t), function(i) 2*pt(abs(result_t$T[i]), result_t$df[i],lower=FALSE)) 
result_t$fdr <- p.adjust(result_t$p_value, "fdr") 
result_t$gene_label <- 1 
result_t$gene_label[result_t$gene == query_gene] <- 2 

result_t_matrix <- result_t %>% dplyr::select(gene, cancer_type, p_value) %>% tidyr::spread(cancer_type, p_value) 

# ordered by difference in means 
head(result_t) 


# ordered by T statistic 
head(result_t[order(result_t$T, decreasing=T), ]) 

# ordered by FDR 
head(result_t[order(result_t$fdr, decreasing=F), ]) 


qplot(data=result_t, x=T, y=mean_diff, shape=as.factor(gene_label), col=as.factor(fdr < 0.05), ylab="Difference in Means", xlab="T statistic") 



