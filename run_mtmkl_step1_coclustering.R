library(AUC)

optimizer <- "mosek"

source("classification_helper.R")
if (optimizer == "cplex") {
  source("solve_classification_models_cplex.R")
}
if (optimizer == "mosek") {
  source("solve_classification_models_mosek.R")
}
source("group_lasso_multitask_multiple_kernel_classification_train.R")
source("group_lasso_multitask_multiple_kernel_classification_test.R")
source("mtmkl_coclustering_algorithm.R")

args <- commandArgs(trailingOnly = TRUE)
replication <- as.numeric(args[[1]]) %% 100 + 1
blocks <- as.numeric(args[[2]]) #number of blocks in a co-clustering

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

full_name <- paste0("TCGA-", paste0(gsub(cohorts, pattern = "TCGA-", replacement = ""), collapse = "-"))
TN <- length(cohorts)

data_path <- "./data"

#The root node of the branching method requires an initial matrix of kernel weights
#We use the kernel weights obtained by Rahimi & Gonen (2018) via gLMKL
initial_data_path <- "initial_data"
initial_info_file <- sprintf("%s/glmkl_hallmark_cohort_weights.RData", initial_data_path)
if(!file.exists(initial_info_file))
{
  stop("Initial info file missing!")
}

branching_iterations <- 20 #maximum depth of the branching tree
branching_time_limit <- 0 #hr*60 minutes, put 0 for unlimited time

result_path <- sprintf("./mtmkl_coclustering_K%d", blocks)

if (dir.exists(sprintf("%s", result_path)) == FALSE) {
  dir.create(sprintf("%s", result_path))
}
if (dir.exists(sprintf("%s/%s", result_path, full_name)) == FALSE) {
  dir.create(sprintf("%s/%s", result_path, full_name)) 
}

state_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_state.RData", result_path, full_name, "hallmark", replication)
result_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_result.RData", result_path, full_name, "hallmark", replication)
prediction_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_prediction.RData", result_path, full_name, "hallmark", replication)
coclustering_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_coclustering.RData", result_path, full_name, "hallmark", replication)

if (file.exists(state_file) == FALSE) {
  epsilon <- 1e-5
  fold_count <- 4
  train_ratio <- 0.8
  iteration_count <- 400 #maximum number of iterations of the cutting-plane (CP) method for each branch
  optimality_gap <- 1e-6 #optimality gap threshold of the CP method performed on each branch
  
  max_svm_cuts <- 500 #Maximum number of initial svm cuts used for warming up the CP method. 
                      #These cust may be obtained and stored in the previous iterations/branches 
  multi_initial_svm_cuts <- FALSE # If TRUE, one cut is added per type. 
                                    #If FALSE, cuts corresponding to diffrent types are integrated into a single cut, to evoid adding too many cuts
  random_seed <- 1606 * replication
  
  coclustering_parameters <- list()
  coclustering_parameters$group_convergance <- 2 #A paramter determing the convergence of the branching process
  coclustering_parameters$branching_iterations <- branching_iterations
  coclustering_parameters$branching_time_limit <- branching_time_limit
  coclustering_parameters$branches <- 2 #number of branches in the heuristic branching process
  coclustering_parameters$diversity_dist <- 2 #a threshold on the diversity among the solutions associated with branches 
  coclustering_parameters$blocks <- blocks
  
  coclustering_parameters$perform_extensive_singel_final_iteration <- TRUE 
        #If TRUE, the best co-clustering solution identified in the branching process 
        #is processed extensively, to obtain exact values of eta, alpha and b
  coclustering_parameters$final_iteration_count <- iteration_count
  coclustering_parameters$final_optimality_gap <- optimality_gap
  
  cutting_plane_parameters <- list()
  cutting_plane_parameters$epsilon <- epsilon
  cutting_plane_parameters$iteration_count <- iteration_count
  cutting_plane_parameters$optimality_gap <- optimality_gap
  cutting_plane_parameters$multi_initial_svm_cuts <- multi_initial_svm_cuts
  cutting_plane_parameters$max_svm_cuts <- max_svm_cuts
  
  X <- vector("list", TN)
  y <- vector("list", TN)
  negative_indices <- vector("list", TN)
  positive_indices <- vector("list", TN)
  for (t in 1:TN) {
    load(sprintf("%s/%s.RData", data_path, cohorts[t]))
    
    common_patients <- intersect(rownames(TCGA$clinical)[which(is.na(TCGA$clinical$pathologic_stage) == FALSE)], rownames(TCGA$mrna))
    
    X[[t]] <- log2(TCGA$mrna[common_patients,] + 1)
    y[[t]] <- rep(NA, length(common_patients))
    
    y[[t]][TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage I",  "Stage IA",  "Stage IB",  "Stage IC")] <- +1
    y[[t]][TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                                     "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                                     "Stage IV",  "Stage IVA",  "Stage IVB",  "Stage IVC")] <- -1
    
    valid_patients <- which(is.na(y[[t]]) == FALSE)
    valid_features <- as.numeric(which(apply(X[[t]][valid_patients,], 2, sd) != 0))
    X[[t]] <- X[[t]][valid_patients, valid_features]
    y[[t]] <- y[[t]][valid_patients]
    
    negative_indices[[t]] <- which(y[[t]] == -1)
    positive_indices[[t]] <- which(y[[t]] == +1)
  }
  
  pathways <- read_pathways("hallmark")
  gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
  for (t in 1:TN) {
    X[[t]] <- X[[t]][, which(colnames(X[[t]]) %in% gene_names)]
  }
  
  P <- length(pathways)
  
  train_negative_indices <- vector("list", TN)
  train_positive_indices <- vector("list", TN)
  negative_allocation <- vector("list", TN)
  positive_allocation <- vector("list", TN)
  for (t in 1:TN) {
    set.seed(random_seed)
    train_negative_indices[[t]] <- sample(negative_indices[[t]], ceiling(train_ratio * length(negative_indices[[t]])))
    train_positive_indices[[t]] <- sample(positive_indices[[t]], ceiling(train_ratio * length(positive_indices[[t]])))
    
    negative_allocation[[t]] <- sample(rep(1:fold_count, ceiling(length(train_negative_indices[[t]]) / fold_count)), length(train_negative_indices[[t]]))
    positive_allocation[[t]] <- sample(rep(1:fold_count, ceiling(length(train_positive_indices[[t]]) / fold_count)), length(train_positive_indices[[t]]))
  }
  
  tuples <- get_cross_validation_tuples_coclustering_step() #returns a list of different combinations of C, lambda_1, and ratio (gamma)
  
  auroc_tuples <- matrix(0, nrow = nrow(tuples), ncol = ncol(tuples)+fold_count) #to be filled with the AUROC values of different folds with different parameters
  colnames(auroc_tuples) <- c("C", "Lambda", "Ratio", paste("Auroc",1:fold_count))
  
  K_train <- vector("list", TN)
  K_test <- vector("list", TN)
  y_train <- vector("list", TN)
  y_test <- vector("list", TN)
  for (fold in 1:fold_count) {
    for (t in 1:TN) {
      train_indices <- c(train_negative_indices[[t]][which(negative_allocation[[t]] != fold)], train_positive_indices[[t]][which(positive_allocation[[t]] != fold)])
      test_indices <- c(train_negative_indices[[t]][which(negative_allocation[[t]] == fold)], train_positive_indices[[t]][which(positive_allocation[[t]] == fold)]) 
      
      X_train <- X[[t]][train_indices,]
      X_test <- X[[t]][test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      
      N_train <- nrow(X_train)
      N_test <- nrow(X_test)
      K_train[[t]] <- array(0, dim = c(N_train, N_train, P))
      K_test[[t]] <- array(0, dim = c(N_test, N_train, P))
      for (m in 1:P) {
        feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
        D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
        D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
        sigma <- mean(D_train)
        K_train[[t]][,,m] <- exp(-D_train^2 / (2 * sigma^2))
        K_test[[t]][,,m] <- exp(-D_test^2 / (2 * sigma^2))
      }
      
      y_train[[t]] <- y[[t]][train_indices]
      y_test[[t]] <- y[[t]][test_indices]
    }
    
    for(tpl in 1:nrow(tuples))
    {
      C <- tuples[tpl, "C"]
      lambda_ratio <- tuples[tpl, "Ratio"]
      lambda_linear <- tuples[tpl, "Lambda"]
      lambda_quadratic <- lambda_linear * lambda_ratio
      
      cutting_plane_parameters$C <- C
      cutting_plane_parameters$lambda_linear <- lambda_linear
      cutting_plane_parameters$lambda_quadratic <- lambda_quadratic
      
      coclustering_parameters$lambda_ratio <- lambda_ratio
      
      var_name <- load(initial_info_file) #loads a file containing the root node eta for branching
      eta_for_branching <- eval(parse(text = var_name))
      
      block_clusters <- generate_block_clusters(eta = eta_for_branching, 
                                                K = blocks, 
                                                lambda_ratio = lambda_ratio, 
                                                min_cluster_size = 2, 
                                                diversity_dist = 2, 
                                                required_diverse_solutions = 4)
      initial_diverse_solutions <- block_clusters$diverse_solutions
      
      outputs <- perform_mtmkl_coclustering_step(K_train, y_train, 
                                                 initial_diverse_solutions, 
                                                 coclustering_parameters, cutting_plane_parameters)
      state <- outputs$state
      prediction <- group_lasso_multitask_multiple_kernel_classification_test(K_test, state$eta, state$alpha, state$b)
      auroc <- numeric(TN)
      for (t in 1:TN) {
        auroc[t] <- auc(roc(prediction$f[[t]], as.factor(y_test[[t]])))
      }
      result <- list(AUROC = auroc)
      auroc_tuples[tpl, c("C", "Ratio", "Lambda", paste("Auroc", fold))] <- c(C, lambda_ratio, lambda_linear, mean(result$AUROC))
    }
  }
  
  average_aurocs <- rowMeans(auroc_tuples[,3+(1:fold_count)])
  tuple_star <- which(average_aurocs == max(average_aurocs))[1]
  C <- auroc_tuples[tuple_star, "C"]
  lambda_ratio <- auroc_tuples[tuple_star, "Ratio"]
  lambda_linear <- auroc_tuples[tuple_star, "Lambda"]
  lambda_quadratic <- lambda_linear * lambda_ratio
  
  for (t in 1:TN) {
    train_indices <- c(train_negative_indices[[t]], train_positive_indices[[t]])
    test_indices <- setdiff(1:length(y[[t]]), train_indices)
    
    X_train <- X[[t]][train_indices,]
    X_test <- X[[t]][test_indices,]
    X_train <- scale(X_train)
    X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
    
    N_train <- nrow(X_train)
    N_test <- nrow(X_test)
    K_train[[t]] <- array(0, dim = c(N_train, N_train, P))
    K_test[[t]] <- array(0, dim = c(N_test, N_train, P))
    for (m in 1:P) {
      feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
      D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
      D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
      sigma <- mean(D_train)
      K_train[[t]][,,m] <- exp(-D_train^2 / (2 * sigma^2))
      K_test[[t]][,,m] <- exp(-D_test^2 / (2 * sigma^2))
    }
    
    y_train[[t]] <- y[[t]][train_indices]
    y_test[[t]] <- y[[t]][test_indices] 
  }
  
  #The tuned parameters are used for the final training:
  
  print(sprintf("running final training: C = %g, lambda_linear = %g, lambda_quadratic = %g", C, lambda_linear, lambda_quadratic))
  
  cutting_plane_parameters$C <- C
  cutting_plane_parameters$lambda_linear <- lambda_linear
  cutting_plane_parameters$lambda_quadratic <- lambda_quadratic
  
  coclustering_parameters$lambda_ratio <- lambda_ratio
  
  var_name <- load(initial_info_file) #loads a file containing the root node eta for branching
  initial_info <- eval(parse(text = var_name))
  eta_for_branching <- initial_info$cohort_weights[,1:TN]
  
  block_clusters <- generate_block_clusters(eta = eta_for_branching, 
                                            K = blocks, 
                                            lambda_ratio = lambda_ratio, 
                                            min_cluster_size = 2, 
                                            diversity_dist = 2, 
                                            required_diverse_solutions = 4)
  initial_diverse_solutions <- block_clusters$diverse_solutions
  
  outputs <- perform_mtmkc_algorithm(K_train, y_train, 
                                     initial_diverse_solutions, 
                                     coclustering_parameters, cutting_plane_parameters)
  
  state <- outputs$state
  state$auroc_tuples <- auroc_tuples
  prediction <- group_lasso_multitask_multiple_kernel_classification_test(K_test, state$eta, state$alpha, state$b)
  auroc <- numeric(TN)
  for (t in 1:TN) {
    auroc[t] <- auc(roc(prediction$f[[t]], as.factor(y_test[[t]])))
  }
  result <- list(AUROC = auroc)
  coclustering_info <- outputs$coclustering_info
  
  save("state", file = state_file)
  save("prediction", file = prediction_file)
  save("result", file = result_file)
  save("coclustering_info", file = coclustering_file)
}