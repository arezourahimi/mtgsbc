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