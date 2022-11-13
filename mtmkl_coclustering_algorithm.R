perform_mtmkl_coclustering_step <- function(K_train, y_train, 
                                            initial_diverse_solutions, 
                                            coclustering_parameters, cutting_plane_parameters)
{
  diverse_solutions <- initial_diverse_solutions
  TN <- length(diverse_solutions[[1]]$cohorts)
  P <- length(diverse_solutions[[1]]$pathways)
  
  start_time <- Sys.time()
  
  final_iteration <- FALSE
  initial_cuts_info <- NULL
  all_initial_cuts <- NULL
  best_branch_groups <- list()
  
  outputs <- list()
  outputs$coclustering_info <- list(tree_info = list(), global_best_branch = list(objective = Inf))
  
  parameters <- list()
  parameters$epsilon <- cutting_plane_parameters$epsilon
  parameters$C <- cutting_plane_parameters$C
  parameters$lambda <- cutting_plane_parameters$lambda_quadratic
  parameters$iteration_count <- cutting_plane_parameters$iteration_count
  parameters$optimality_gap <- cutting_plane_parameters$optimality_gap
  
  iteration <- 1
  
  while(TRUE)
  {
    current_iteration_branches <- coclustering_parameters$branches
    if(final_iteration & coclustering_parameters$perform_extensive_singel_final_iteration)
    {
      parameters$iteration_count <- coclustering_parameters$final_iteration_count
      parameters$optimality_gap <- coclustering_parameters$final_optimality_gap
      current_iteration_branches <- 1
    }
    
    best_branch <- list(objective = Inf)
    outputs$coclustering_info$tree_info[[iteration]] <- list()
    
    for(d in 1:current_iteration_branches)
    {
      print(sprintf("CC: Iteration %d: blocks = %d, branch = %d (of %d) distance = %d", iteration, coclustering_parameters$blocks, d, 
                    current_iteration_branches, coclustering_parameters$diversity_dist))
      
      block_penalty_info <- generate_block_penalty_info_from_coclusters(diverse_solutions[[d]])
      
      penalty_vector <- block_penalty_info$penalty_vector * cutting_plane_parameters$lambda_linear
      penalty_matrix <- block_penalty_info$penalty_matrix
      
      output <- group_lasso_multitask_multiple_kernel_classification_cutting_plane_train(K_train, y_train, parameters, 
                                                                                         penalty_matrix, penalty_vector,
                                                                                         initial_cuts_info)
      state <- output$state
      
      if(!is.null(output$new_cuts_info))
      {
        if(!cutting_plane_parameters$multi_initial_svm_cuts)
        {
          cuts <- length(output$new_cuts_info$SVM_cuts$rhs) / TN
          new_svm_constraints <- matrix(0, nrow = cuts, ncol = ncol(output$new_cuts_info$SVM_cuts$matrix))
          new_svm_rhs <- numeric(cuts)
          
          for(cut in 1:cuts)
          {
            cut_range = (1:TN) + TN*(cut-1)
            new_svm_constraints[cut,] <- colSums(output$new_cuts_info$SVM_cuts$matrix[cut_range,])
            new_svm_rhs[cut] <- sum(output$new_cuts_info$SVM_cuts$rhs[cut_range])
          }
        }else
        {
          new_svm_constraints <- output$new_cuts_info$SVM_cuts$matrix
          new_svm_rhs <- output$new_cuts_info$SVM_cuts$rhs
        }
        
        new_etas <- output$new_cuts_info$etas
        
        if(is.null(all_initial_cuts))
        {
          all_initial_cuts <- list(SVM_cuts = list(matrix = new_svm_constraints, rhs = new_svm_rhs))
        }else
        {
          all_initial_cuts$SVM_cuts$matrix = rbind(all_initial_cuts$SVM_cuts$matrix, new_svm_constraints)
          all_initial_cuts$SVM_cuts$rhs = c(all_initial_cuts$SVM_cuts$rhs, new_svm_rhs)
        }
        
        initial_cuts_info = all_initial_cuts
        
        sc <- length(initial_cuts_info$SVM_cuts$rhs)
        if(sc > max_svm_cuts)
        {
          initial_cuts_info$SVM_cuts$rhs <- initial_cuts_info$SVM_cuts$rhs[(sc-max_svm_cuts+1):sc]
          initial_cuts_info$SVM_cuts$matrix <- initial_cuts_info$SVM_cuts$matrix[(sc-max_svm_cuts+1):sc,]
        }
      }
      
      branch_info <- list()
      branch_info$objective <- state$iteration_info[nrow(state$iteration_info),2]
      branch_info$alpha <- state$alpha
      branch_info$b <- state$b
      
      eta_matrix <- matrix(0, nrow = P, ncol = TN)
      colnames(eta_matrix) <- cohorts_without_TCGA
      for (t in 1:TN) {
        eta_matrix[,t] <- state$eta[[t]]
      }
      branch_info$eta <- eta_matrix
      branch_info$block_clusters <- diverse_solutions[[d]]
      
      deviation_cohort <- Inf
      deviation_pathway <- Inf
      if(length(best_branch_groups) > 0)
      {
        for(g in best_branch_groups)
        {
          dc <- calculate_cluster_distance(g$cohorts, diverse_solutions[[d]]$cohorts, blocks)
          dp <- calculate_cluster_distance(g$pathways, diverse_solutions[[d]]$pathways, blocks)
          if(dc < deviation_cohort)
          {
            deviation_cohort <- dc
            deviation_pathway <- dp
          }else if(dc == deviation_cohort)
          {
            if(dp < deviation_pathway)
              deviation_pathway <- dp
          }
        }
      }else
      {
        deviation_cohort <- TN
        deviation_pathway <- P
      }
      branch_info$deviation_cohort <- deviation_cohort
      branch_info$deviation_pathway <- deviation_pathway
      
      outputs$coclustering_info$tree_info[[iteration]][[d]] <- list(objective = branch_info$objective,
                                                                    block_clusters = branch_info$block_clusters)
      
      if(branch_info$objective < best_branch$objective)
      {
        best_branch <- branch_info
        
        if(branch_info$objective < outputs$coclustering_info$global_best_branch$objective)
        {
          outputs$coclustering_info$global_best_branch <- branch_info
          outputs$state <- state
        }
      }
      
      print(sprintf("*** CC: Iteration %d branch %d done with objective = %f, grouping deviation: c=%d, p=%d", 
                    iteration, d, branch_info$objective, deviation_cohort, deviation_pathway))
      
      print("cohort groups:")
      for(k in 1:blocks) 
        print(sprintf(" %d: %s",k, paste0(cohorts_without_TCGA[which(diverse_solutions[[d]]$cohorts ==k)], collapse = ", ")))
      
      print("pathway groups:")
      for(k in 1:blocks) 
        print(sprintf(" %d: %s",k, paste0(c(1:P)[which(diverse_solutions[[d]]$pathways ==k)], collapse = ", ")))
    }
    
    if(!is.null(best_branch$eta))
    {
      if(final_iteration)
      {
        break
      }
      
      deviation <- best_branch$deviation_cohort + best_branch$deviation_pathway
      
      if(deviation <= coclustering_parameters$group_convergance)
      {
        final_iteration <- TRUE
      }
      
      if(coclustering_parameters$branching_iterations > 0 & iteration >= coclustering_parameters$branching_iterations)
      {
        final_iteration <- TRUE
      }
      
      if(coclustering_parameters$branching_time_limit > 0)
      {
        total_elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        if(total_elapsed_time >= branching_time_limit)
        {
          print(sprintf("Time elapsed: %f, terminated at iteration %d due to time limit of %d", total_elapsed_time, iteration, coclustering_parameters$branching_time_limit))
          final_iteration <- TRUE
        }
      }
      
      if(final_iteration & !coclustering_parameters$perform_extensive_singel_final_iteration)
      {
        break
      }
      
      best_branch_groups[[length(best_branch_groups) + 1]] = best_branch$block_clusters
      
      eta_for_branching <- best_branch$eta
      block_clusters <- generate_block_clusters(eta = eta_for_branching, 
                                                K = coclustering_parameters$blocks, 
                                                lambda_ratio = coclustering_parameters$lambda_ratio, 
                                                min_cluster_size = 2, 
                                                diversity_dist = coclustering_parameters$diversity_dist, 
                                                required_diverse_solutions = 2*coclustering_parameters$branches)
      
      diverse_solutions <- block_clusters$diverse_solutions
      
      iteration <- iteration + 1
    }else
    {
      break
    }
  }
  
  return(outputs)
}

#Given a matrix of kernel weights, the number of clusters along with the feasibility requirements,
#returns a near optimal coclustering solution. To prevent stagnation at local optima, 
#it is used as a sub-module for the global genetic algorithm. 
coclustering_heuristic_algorithm_local <- function(eta, K, lambda_1, lambda_2, 
                                                 min_group_size = 2, max_group_size_ratio = 1,
                                                 mutations = 3, swaps = 1, max_nonimproving = 10000,
                                                 initial_solution = NULL, fix_cohorts = FALSE, 
                                                 show_results = FALSE){
  
  TN <- ncol(eta)
  P <- nrow(eta)
  
  Delta = matrix(0, nrow = TN, ncol = TN)
  
  for(s in 1:(TN-1))
  {
    for(q in (s+1):TN)
    {
      Delta[s,q] = Delta[q,s] = sum((eta[,s]-eta[,q])^2)
    }
  }
  
  z <- function(TC, PC){
    obj <- 0
    if(lambda_1 > 0)
    {
      obj1 <- 0
      for(m in 1:P)
      {
        for(s in 1:TN)
        {
          if(TC[s] != PC[m])
          {
            obj1 <- obj1 + eta[m,s]
          }
        }
      }
      
      obj <- obj + lambda_1*obj1
    }
    
    if(lambda_2 > 0)
    {
      obj2 <- 0
      for(s in 1:(TN-1)){
        for(q in (s+1):TN){
          if(TC[s] == TC[q]){
            obj2 <- obj2 + Delta[s,q]
          }
        }
      }
      
      obj <- obj + lambda_2*2*obj2
    }
    
    return(obj)
  }
  
  get_size <- function(group_lables, k)
  {
    return(length(which(group_lables==k)))
  }
  
  group_size_range <- function(group_labels)
  {
    min_size <- length(group_labels)
    max_size <- 0
    for(k in 1:K)
    {
      size_k <- get_size(group_labels, k)
      if(size_k < min_size)
        min_size = size_k
      if(size_k > max_size)
        max_size = size_k
    }
    
    return(c(min_size, max_size))
  } 
  
  if(is.null(initial_solution))
  {
    # random initial sol:
    while(TRUE)
    {
      TC <- sample(1:K,TN, replace = TRUE)
      sizes <- group_size_range(TC)
      if(sizes[1] >= min_group_size & sizes[2] <= max_group_size_ratio*TN)
        break
    }
    
    while(TRUE)
    {
      PC <- sample(1:K, P, replace = TRUE)
      sizes <- group_size_range(PC)
      if(sizes[1] >= min_group_size & sizes[2] <= max_group_size_ratio*P)
        break
    }
    
    obj <- z(TC,PC)
    
    incumbent <- list(cohorts = TC, pathways = PC, objective = obj)
  }else
  {
    incumbent <- initial_solution
  }
  
  iteration <- 1
  nonimproving <- 0
  
  incumbent_solutions = list()
  incumbent_solutions[[1]] = incumbent
  
  while(TRUE){
    
    TC_ <- incumbent[[1]]
    PC_ <- incumbent[[2]]
    Z_ <- incumbent[[3]]
    
    # mutate
    if(mutations > 0)
    {
      for(j in 1:mutations){
        if(!fix_cohorts)
        {
          while(TRUE)
          {
            s <- sample(1:TN,1)
            TC_[s] <- sample(1:K,1)
            sizes <- group_size_range(TC_)
            if(sizes[1] >= min_group_size & sizes[2] <= max_group_size_ratio*TN)
              break
          }
        }
        
        while(TRUE)
        {
          m <- sample(1:P,1)
          PC_[m] <- sample(1:K,1)
          sizes <- group_size_range(PC_)
          if(sizes[1] >= min_group_size & sizes[2] <= max_group_size_ratio*P)
            break
        }
      }
    }
    
    Z_ <- z(TC_,PC_)
    
    #swap
    if(swaps>0)
    {
      for(j in 1:swaps)
      {
        if(!fix_cohorts)
        {
          s1 <- sample(1:TN,1)
          sk1 <- TC_[s1]
          s2 <- sample(1:TN,1)
          sk2 <- TC_[s2]
          while(sk1 == sk2)
          {
            s2 <- sample(1:TN,1)
            sk2 <- TC_[s2]
          }
          
          TC_[s1] <- sk2
          TC_[s2] <- sk1
        }
        
        m1 <- sample(1:P,1)
        mk1 <- PC_[m1]
        m2 <- sample(1:P,1)
        mk2 <- PC_[m2]
        while(mk1 == mk2)
        {
          m2 <- sample(1:P,1)
          mk2 <- PC_[m2]
        }
        
        PC_[m1] <- mk2
        PC_[m2] <- mk1
        
        Z__ <- z(TC_,PC_)
        if(Z__ < Z_)
        {
          Z_ <- Z__
        }else
        {
          if(!fix_cohorts)
          {
            TC_[s1] <- sk1
            TC_[s2] <- sk2
          }
          PC_[m1] <- mk1
          PC_[m2] <- mk2
        }
      }
    }
    
    if(Z_ < incumbent[[3]]){
      incumbent <- list(cohorts = TC_, pathways = PC_, objective = Z_)
      nonimproving <- 0
      
      incumbent_solutions[[length(incumbent_solutions)+1]] <- incumbent
    }else
    {
      nonimproving <- nonimproving + 1
    }
    
    if(show_results)
    {
      if(iteration %% 10000 == 0)
        print(sprintf("iteration %d: obj: %f, incumbent: %f", iteration, Z_, incumbent[[3]]))
    }
    
    iteration <- iteration + 1
    if(nonimproving > max_nonimproving)
    {
      if(show_results)
        print(sprintf("iteration %d: obj: %f, incumbent: %f", iteration, Z_, incumbent[[3]]))
      break
    }
  }
  
  return(incumbent_solutions)
}

#The genetic algorithm is restarted multiple times to decrease the chance premature convergence 
coclustering_heuristic_algorithm <- function(eta, K, lambda_1, lambda_2, 
                                                  min_group_size, max_group_size_ratio,
                                                  mutations, swaps, 
                                                  max_local_nonimproving = 10000, 
                                                  max_global_nonimproving = 10, 
                                                  diversity_dist = 10, required_diverse_solutions = 10, reoptimize = TRUE)
{
  incumbent = list(cohorts = NULL, pathways = NULL, objective = Inf)
  nonimproving <- 0
  replication <- 1
  all_incumbent_solutions = c()
  
  while(TRUE)
  {
    rep_solutions <- coclustering_heuristic_algorithm_local(eta = eta, K = K, 
                                                          lambda_1 = lambda_1, lambda_2 = lambda_2, 
                                                          min_group_size = min_group_size, max_group_size_ratio = max_group_size_ratio,
                                                          mutations = mutations, swaps = swaps, 
                                                          max_nonimproving = max_local_nonimproving, show_results = FALSE, initial_solution = NULL, fix_cohorts = FALSE)
    
    all_incumbent_solutions = c(all_incumbent_solutions, rep_solutions)
    rep <- rep_solutions[[length(rep_solutions)]]
    if(rep[[3]] < incumbent[[3]])
    {
      incumbent <- rep
      nonimproving <- 0
    }else
    {
      nonimproving <- nonimproving + 1
    }
    
    print(sprintf("replication %d: obj: %f, incumbent: %f", replication, rep[[3]], incumbent[[3]]))
    
    if(nonimproving >= max_global_nonimproving)
    {
      break
    }
    replication <- replication + 1
  }
  
  diverse_solutions <- select_diverse_clusters(all_incumbent_solutions, diversity_dist, required_diverse_solutions)
  
  if(reoptimize)
  {
    print("reoptimizing diverse solutions:")
    for(d in 1:length(diverse_solutions))
    {
      ds <- diverse_solutions[[d]]
      opt <- coclustering_heuristic_algorithm_local(eta = eta, K = K, 
                                                  lambda_1 = lambda_1, lambda_2 = lambda_2, 
                                                  min_group_size = min_group_size, max_group_size_ratio = max_group_size_ratio, 
                                                  mutations = mutations, swaps = swaps, 
                                                  max_nonimproving = max_local_nonimproving, fix_cohorts = TRUE, initial_solution = ds, 
                                                  show_results = FALSE)
      ds_opt <- opt[[length(opt)]]
      print(sprintf("%d: objective reduced from %f to %f",d , ds$objective, ds_opt$objective))
      diverse_solutions[[d]] <- ds_opt
    }
    
    diverse_solutions <- diverse_solutions[order(sapply(diverse_solutions, function(incumbent) incumbent$objective))]
  }
  
  return(diverse_solutions)
}

#Given a list of solutions, returns the top diverse solutions. A minimum diversity is enforced on each pair of the selected solutions
select_diverse_clusters <- function(all_incumbent_solutions, diversity_dist, required_diverse_solutions)
{
  sorted_incumbent_solutions <- all_incumbent_solutions[order(sapply(all_incumbent_solutions, function(incumbent) incumbent$objective))]
  
  diverse_solutions <- list()
  K <- length(unique(sorted_incumbent_solutions[[1]]$cohorts))
  
  for(inc in sorted_incumbent_solutions)
  {
    if(length(diverse_solutions) == 0)
    {
      diverse_solutions[[1]] <- inc
    }else
    {
      min_dist <- Inf
      for(d in 1:length(diverse_solutions))
      {
        dist <- calculate_cluster_distance(inc$cohorts, diverse_solutions[[d]]$cohorts, K)
        if(dist < min_dist)
          min_dist <- dist
      }
      
      if(min_dist >= diversity_dist)
      {
        index <- length(diverse_solutions) + 1
        diverse_solutions[[index]] <- inc
      }
    }
    
    if(length(diverse_solutions) >= required_diverse_solutions)
      break()
  }
  
  return(diverse_solutions)
}