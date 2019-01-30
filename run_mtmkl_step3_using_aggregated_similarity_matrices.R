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

args <- commandArgs(trailingOnly = TRUE)
replication <- as.numeric(args[[1]]) %% 100 + 1

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

full_name <- paste0("TCGA-", paste0(gsub(cohorts, pattern = "TCGA-", replacement = ""), collapse = "-"))
TN <- length(cohorts)

data_path <- "./data"
initial_data_path <- "initial_data" #this is the directory to a file containing the similarity matrices produced as a result of step 2

result_path <- sprintf("./result_mtmkl")

if (dir.exists(sprintf("%s", result_path)) == FALSE) {
  dir.create(sprintf("%s", result_path)) 
}
if (dir.exists(sprintf("%s/%s", result_path, full_name)) == FALSE) {
  dir.create(sprintf("%s/%s", result_path, full_name)) 
}

#Before perfming the final step, we need the similarity matrices that were obtained by aggregating the results of the coclustering step
similarity_matrices_file <- sprintf("%s/similarity_matrices_hallmark.RData", initial_data_path)
if(!file.exists(similarity_matrices_file))
{
  stop("Similarity matrices are missing! Use run_mtmkl_step2_collect_coclustering_data.R for generating the similarity matrices.")
}

state_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_state.RData", result_path, full_name, "hallamrk", replication)
result_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_result.RData", result_path, full_name, "hallamrk", replication)
prediction_file <- sprintf("%s/%s/glmtmklp_pathway_%s_measure_AUROC_replication_%d_prediction.RData", result_path, full_name, "hallamrk", replication)

if (file.exists(state_file) == FALSE) {
  epsilon <- 1e-5
  fold_count <- 4
  train_ratio <- 0.8
  iteration_count <- 400
  optimality_gap <- 1e-6
  random_seed <- 1505 * replication
  
  var_name <- load(similarity_matrices_file)
  similarity_matrices <- eval(parse(text = var_name))
  
  penalty_info <- generate_block_penalty_info(similarity_matrices)
  penalty_vector <- penalty_info$penalty_vector
  penalty_matrix <- penalty_info$penalty_matrix
  
  parameters <- list()
  parameters$epsilon <- epsilon
  parameters$iteration_count <- iteration_count
  parameters$optimality_gap <- optimality_gap
  
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
  
  pathways <- read_pathways("hallamrk")
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
  
  tuples <- get_cross_validation_tuples_mtmkl()
  
  auroc_tuples <- matrix(0, nrow = nrow(tuples), ncol = 3+fold_count)
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
      
      print(sprintf("running fold = %d, C = %g, lambda_linear = %g, lambda_quadratic = %g", fold, C, lambda_linear, lambda_quadratic))
      
      parameters$C <- C
      parameters$lambda <- lambda_quadratic
      
      pv <- penalty_vector * lambda_linear
      output <- group_lasso_multitask_multiple_kernel_classification_cutting_plane_train(K_train, y_train, 
                                                                                         parameters, 
                                                                                         penalty_matrix, penalty_vector = pv)
      state <- output$state
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
  
  print(sprintf("running final training: C = %g, lambda_linear = %g, lambda_quadratic = %g", C, lambda_linear, lambda_quadratic))
  
  parameters$iteration_count <- iteration_count
  parameters$optimality_gap <- optimality_gap
  parameters$C <- C
  parameters$lambda <- lambda_quadratic
  
  pv <- penalty_vector * lambda_linear
  output <- group_lasso_multitask_multiple_kernel_classification_cutting_plane_train(K_train, y_train, 
                                                                                     parameters, 
                                                                                     penalty_matrix, penalty_vector = pv)
  state <- output$state
  state$auroc_tuples <- auroc_tuples
  prediction <- group_lasso_multitask_multiple_kernel_classification_test(K_test, state$eta, state$alpha, state$b)
  auroc <- numeric(TN)
  for (t in 1:TN) {
    auroc[t] <- auc(roc(prediction$f[[t]], as.factor(y_test[[t]])))
  }
  result <- list(AUROC = auroc)
  
  save("state", file = state_file)
  save("prediction", file = prediction_file)
  save("result", file = result_file)
}
