source("classification_helper.R")

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

cohorts_without_TCGA <- gsub(cohorts, pattern = "TCGA-", replacement = "")
full_name <- paste0("TCGA-", paste0(cohorts_without_TCGA, collapse = "-"))
TN <- length(cohorts)

result_path <- sprintf("./result_mtmkl")

significance <- 0.01 #used for distinguishing the significant kernel weights across different types and pathways

pathways <- read_pathways("hallmark")
P <- length(pathways)

pathway_names <- sapply(pathways, function(p) {p$name})
if(startsWith(tolower(pathway_names[1]), "hallmark"))
  pathway_names <- substr(pathway_names, nchar("hallmark") + 2, nchar(pathway_names))

replication_count <- 100

cohort_weights <- matrix(0, P, TN) #will be populated by the average kerenl weights of different cohorts
colnames(cohort_weights) <- cohorts_without_TCGA
rownames(cohort_weights) <- pathway_names

cohort_significance <- matrix(0, P, TN) #will be populated by the relative frequency of significant kerenl weights of different cohorts
colnames(cohort_significance) <- cohorts_without_TCGA
rownames(cohort_significance) <- pathway_names

aurocs <- matrix(0, nrow = replication_count, ncol = TN)
colnames(aurocs) <- cohorts_without_TCGA

loaded_replications <- 0

for (replication in 1:replication_count) {
  
  state_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_state.RData", result_path, full_name, "hallmark", replication)
  result_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_result.RData", result_path, full_name, "hallmark", replication)
  
  if(file.exists(state_file))
  {
    loaded_replications <- loaded_replications + 1
    
    load(state_file)
    load(result_file)
    
    etas = state$eta
    
    names(etas) = cohorts_without_TCGA
    
    aurocs[replication, ] = result$AUROC
    
    for (cohort in cohorts_without_TCGA) {
      
      eta  = etas[[cohort]]
      
      cohort_weights[,cohort] <- cohort_weights[,cohort] + eta
      cohort_significance[,cohort] <- cohort_significance[,cohort] + 1*(eta>significance)
    }
  }
}

cohort_weights <- cohort_weights/loaded_replications 
cohort_significance <- cohort_significance/loaded_replications 

mtmkl_info <- list(cohort_weights = cohort_weights, 
                   cohort_significance = cohort_significance, 
                   significance = significance,
                   aurocs = aurocs)

mtmkl_info_file <- sprintf("%s/mtmkl_combined_results.RData", result_path)

save("mtmkl_info", file = mtmkl_info_file)