#execute this code to collect/aggregate the outpus of the co-clusterinng step performed over 
#the 100 replications, before executing the finanl step of the MTMKL algorithm

source("classification_helper.R")

args <- commandArgs(trailingOnly = TRUE)
blocks <- as.numeric(args[[1]])

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

cohorts_without_TCGA <- gsub(cohorts, pattern = "TCGA-", replacement = "")
TN <- length(cohorts)
full_name <- paste0("TCGA-", paste0(cohorts_without_TCGA, collapse = "-"))

result_path <- sprintf("./mtmkl_coclustering_K%d", blocks)

pathways <- read_pathways("hallmark")
P <- length(pathways)
pathway_names <- sapply(pathways, function(p) {p$name})
if(startsWith(tolower(pathway_names[1]), "hallmark"))
  pathway_names <- substr(pathway_names, nchar("hallmark") + 2, nchar(pathway_names))

replication_info <- list()

cohort_cohort_frequency <- matrix(0, nrow = TN, ncol = TN)
colnames(cohort_cohort_frequency) <- cohorts_without_TCGA
rownames(cohort_cohort_frequency) <- cohorts_without_TCGA

pathway_pathway_frequency <- matrix(0, nrow = P, ncol = P)
rownames(pathway_pathway_frequency) <- pathway_names
colnames(pathway_pathway_frequency) <- pathway_names

pathway_cohort_frequency = matrix(0, nrow = P, ncol = TN)
colnames(pathway_cohort_frequency) <- cohorts_without_TCGA
rownames(pathway_cohort_frequency) <- pathway_names

replication_count <- 100
loaded_replications <- 0

for(replication in 1:replication_count)
{
  coclustering_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_coclustering.RData", result_path, full_name, "hallmark", replication)
  
  if(file.exists(coclustering_file))
  {
    load(coclustering_file)
    loaded_replications <- loaded_replications + 1
    
    final_clustering <- coclustering_info$global_best_branch$block_clusters
    
    cg <- final_clustering$cohorts
    pg <- final_clustering$pathways
    cohort_cohort_frequency <- cohort_cohort_frequency + calculate_cocluster_matrix(cg)
    pathway_pathway_frequency <- pathway_pathway_frequency + calculate_cocluster_matrix(pg)
    pathway_cohort_frequency <- pathway_cohort_frequency + calculate_cocluster_matrix(pg, cg)
  }
}

if(loaded_replications == 0)
{
  stop("Failed to locate the outputs of coclustering step!")
}

cohort_cohort_frequency <- cohort_cohort_frequency / loaded_replications
pathway_pathway_frequency <- pathway_pathway_frequency / loaded_replications
pathway_cohort_frequency <- pathway_cohort_frequency / loaded_replications

similarity_matrices = list(cohorts = cohort_cohort_frequency, 
                                      pathways = pathway_pathway_frequency, 
                                      pathway_cohorts = pathway_cohort_frequency)

initial_data_path <- "initial_data"
similarity_matrices_file <- sprintf("%s/similarity_matrices_hallmark.RData", initial_data_path)
save(similarity_matrices, file = similarity_matrices_file)