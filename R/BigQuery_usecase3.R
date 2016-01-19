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
#' In this query, we get that a fair number of Samples in BRCA are mutated for PIK3CA which make this gene
#' a suitable candidate for our test
#' 
#' Now we are going to compare the expression for all genes in the mutated samples vs the samples that don't
#' have a mutation for the query gene. We are going to group the data by gene and study (cancer type)
#'  
## ----chunk5-------------------------------------------------------------- 
ptm1 <- proc.time() 


# Now we are ready to run the t-test query query. 
result = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "t_test_viral_status_expression.sql"), 
  project=project #, 
#  replacements=list() 
  
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
