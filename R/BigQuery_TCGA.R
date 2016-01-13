#'
#'source("http://bioconductor.org/biocLite.R")
#'
#'
#' # Exploring the TCGA data in BigQuery 
#'  
#' The ISB-CGC (isb-cgc.org) project has aggregated and curated all of the TCGA open-access clinical, biospecimen, and Level-3 molecular data and uploaded it into BigQuery tables that are open to the public.  Here we will show you how you can begin to work with these tables from the familiar R environment. 
#'  
#' ### Helpful BigQuery links 
#'  
#' For this example, we'll also be working with [Google BigQuery](https://cloud.google.com/bigquery/). It's often helpful to have a [link to the docs](https://cloud.google.com/bigquery/what-is-bigquery) handy, and especially the [query reference](https://cloud.google.com/bigquery/query-reference). 
#'  
#' ## Run a query from R 
#'  
#' We will start by loading four R packages: 
#' - the [bigrquery](https://github.com/hadley/bigrquery) package written by Hadley Wickham implements an R interface to [Google BigQuery](https://cloud.google.com/bigquery/), 
#' - the [dplyr](https://github.com/hadley/dplyr) package provides a set of tools for efficiently manipulating datasets in R, and 
#' - the [ggplot2](https://github.com/hadley/ggplot2) package for elegant graphics, and 
#' - the [scales](https://github.com/hadley/scales) package for visualization-oriented scale functions. 
#'  
## ----message=FALSE------------------------------------------------------- 
require(dplyr) || install.packages("dplyr") 
require(bigrquery) || install.packages("bigrquery") 
require(scales) || install.packages("scales") 
require(ggplot2) || install.packages("ggplot2") 
#library(ISBCGCExamples)

# The directory in which the files containing SQL reside.
sqlDir = file.path("Path to the sql files")
#'  
## ----eval=FALSE---------------------------------------------------------- 
## ######################[ TIP ]######################################## 
## ## Set the Google Cloud Platform project id under which these queries will run. 
## ## 
## ## If you are using the Google Bioconductor workshop docker image, this is already 
## ## set for you in your .Rprofile and you can skip this step. 
##  
## # project <- "YOUR-PROJECT-ID" 
## ##################################################################### 
 
#'  
#' ## Spearman Correlation in BigQuery 
#'  
#' We will start by performing the correlation directly in BigQuery.  We will use a pre-defined SQL query in which key strings will be replaced according to the values specified in the "replacements" list. 
#'  
## ----comment=NA---------------------------------------------------------- 
# Set the desired tables to query. 
expressionTable <- "isb-cgc:tcga_201510_alpha.mRNA_UNC_HiSeq_RSEM" 
query_gene <- "EGFR"
# Do not correlate unless there are at least this many observations available 
minimumNumberOfObservations = 30 


#'  
#' Note that when you send the first query, you will need to go through the authentication flow with BigQuery.  You will be provided with a url to cut and  paste into your browser, and then you will get an authorization code to cut and paste back here. 
#'  
#' [bigrquery](https://github.com/hadley/bigrquery) uses the package [httr](https://github.com/hadley/httr) to perform OAuth. 
#'  
#' Let's start by just counting the number of records in the table. 
#' First we'll just create the query string and echo it: 
## ------------------------------------------------------------------------ 
querySql <- paste("SELECT COUNT(1) as total_rows FROM [", expressionTable, "]", sep="") 
querySql 

#'  
#' And then we'll send the query to the cloud for execution: 
## ------------------------------------------------------------------------ 
result <- query_exec(querySql, project=project) 
result 

  
#' We'll pause for a moment here and have a look at the size of our cohort table: 
#'  
## ----comment=NA---------------------------------------------------------- 
cohortInfo <- get_table("isb-cgc","tcga_201510_alpha","mRNA_UNC_HiSeq_RSEM") 
cohortInfo$description 
cohortInfo$numRows 
 
ptm1 <- proc.time() 
 

# Now we are ready to run the spearman correlation query. 
result = DisplayAndDispatchQuery ( 
              file.path(sqlDir, "PearsonCorrelation.sql"), 
              project=project, 
              replacements=list("_EXPRESSION_TABLE_"=expressionTable, 
                                "_QUERY_GENE_"= query_gene, 
                                "_MINIMUM_NUMBER_OF_OBSERVATIONS_"=minimumNumberOfObservations)
                                ) 


ptm2 <- proc.time() - ptm1 
cat("Wall-clock time for BigQuery:",ptm2[3]) 
 
#' Number of rows returned by this query: `nrow(result)`. 
#'  
#' The result is a table with one row for each (gene,protein) pair for which at least 30 data values exist for the specified cohort.  The (gene,protein) pair is defined by a gene symbol and a protein name.  In many cases the gene symbol and the protein name may be identical, but for some genes the RPPA dataset may contain expression values for more than one post-translationally-modified protein product from a particular gene. 
#'  
## ------------------------------------------------------------------------ 
head(result) 


#'  
## ----pearson_density, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA---- 
library(ggplot2) 

# Use ggplot to create a histogram overlaid with a transparent kernel density curve 
ggplot(result, aes(x=pearson_corr)) + 
      # use 'density' instead of 'count' for the histogram 
      geom_histogram(aes(y=..density..), 
                     binwidth=.05, 
                     colour="black", fill="white") + 
      # and overlay with a transparent density plot 
      geom_density(alpha=.2, fill="#FF6666") + 
      # and add a vertical line at x=0 to emphasize that most correlations are positive 
      geom_vline(xintercept=0) 



#'  
#' The R-based Spearman correlation matches the BigQuery result. 
#'  
#' ## Provenance 
#' 
#' 
#' 

ptm3 <- proc.time() 

# Now we are ready to run the pearson correlation query. 
result_spearman = DisplayAndDispatchQuery ( 
  file.path(sqlDir, "SpearmanCorrelation.sql"), 
  project=project, 
  replacements=list("_EXPRESSION_TABLE_"=expressionTable, 
                    "_QUERY_GENE_"= query_gene, 
                    "_MINIMUM_NUMBER_OF_OBSERVATIONS_"=minimumNumberOfObservations)
) 

ptm4 <- proc.time() - ptm3 
cat("Wall-clock time for BigQuery:",ptm4[3]) 


head(result_spearman)
#'  
## ----spearman_density, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA---- 
library(ggplot2) 

# Use ggplot to create a histogram overlaid with a transparent kernel density curve 
ggplot(result_spearman, aes(x=spearman_corr)) + 
  # use 'density' instead of 'count' for the histogram 
  geom_histogram(aes(y=..density..), 
                 binwidth=.05, 
                 colour="black", fill="white") + 
  # and overlay with a transparent density plot 
  geom_density(alpha=.2, fill="#FF6666") + 
  # and add a vertical line at x=0 to emphasize that most correlations are positive 
  geom_vline(xintercept=0) 



## ----provenance, comment=NA---------------------------------------------- 
sessionInfo() 




