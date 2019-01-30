group_lasso_multitask_multiple_kernel_classification_cutting_plane_train <- function(Km, y, parameters,
                                                                                     penalty_matrix, penalty_vector = NULL, 
                                                                                     initial_cuts_info = NULL) {
  # parameters should contain: lambda (quadratic), optimality_gap
  
  TN <- length(Km)
  P <- dim(Km[[1]])[3]
  
  terminate <- FALSE
  UB <- Inf
  LB <- -Inf
  
  #incumbent solution
  incumbent_solution <- list(eta = NULL, alphs = NULL, b = NULL)
  
  iteration <- 1
  
  #define constraint matrix here, and add the equality constraint. Also define the objective function
  #order of variables: Gamma_1, ..., Gamma_T, eta_1_1,...,eta_1_P,eta_2_1,...,eta_2_P,...,eta_T_1,...,eta_T_P
  NVars <- TN + TN*P
  master_obj_coef <- numeric(NVars)
  master_obj_coef[1:TN] <- 1
  if(!is.null(penalty_vector))
  {
    master_obj_coef[(TN+1):NVars] <- penalty_vector
  }
  
  if(!is.null(penalty_matrix))
  {
    obj_quadratic <- matrix(0, NVars, NVars)
    obj_quadratic[(TN+1):NVars, (TN+1):NVars] <- parameters$lambda*penalty_matrix 
  }else
  {
    obj_quadratic <- NULL
  }
  
  master_constraint_matrix <- matrix(0, nrow = TN, ncol = NVars)
  for(t in 1:TN)
  {
    master_constraint_matrix[t,(TN+1+(t-1)*P):(TN+1+(t)*P-1)] <- 1
  }
  
  master_rhs <- rep(1,TN)
  master_senses <- rep("E",TN)
  master_lb <- numeric(NVars)
  
  if(!is.null(initial_cuts_info))
  {
    #these initial cuts may be provided from the former replicaitons of the algorithm
    if(!is.null(initial_cuts_info$SVM_cuts))
    {
      master_constraint_matrix <- rbind(master_constraint_matrix, initial_cuts_info$SVM_cuts$matrix)
      master_rhs <- c(master_rhs, initial_cuts_info$SVM_cuts$rhs)
      master_senses <- c(master_senses, rep("G", length(initial_cuts_info$SVM_cuts$rhs)))
      
      print(sprintf("Added %d initial SVM cuts.", length(initial_cuts_info$SVM_cuts$rhs)))
    }
  }

  #we store the newly generated cuts, in case they are needed in the future  
  new_svm_constraints <- matrix(0, nrow = 0, ncol = NVars)
  new_svm_rhs <- c()
  
  J <- numeric(TN)
  Gamma <- numeric(TN) #approximator of J
  
  alpha <- vector("list", TN)
  b <- vector("list", TN)
  eta <- vector("list", TN)
  for (t in 1:TN) {
    eta[[t]] <- rep(1 / P, P)
  }
  
  quadratic_penalty <- 0 #the (quadratic) penalty as calculated in the master problem
  
  objectives <- c()
  
  info_cols <- c("obj", "UB", "LB", "Gap", "Gamma", "J", "Q-Penalty", "L-Penalty")
  iteration_info <- matrix(0,0,length(info_cols))
  colnames(iteration_info) <- info_cols
  
  while(terminate == FALSE)
  {
    #SOLVE MP and obtain Gamma[t], quadratic_penalty, and eta[[t]]
    master_result <- solve_multitask_quadratic_master_problem(master_obj_coef, obj_quadratic, master_constraint_matrix, master_rhs, master_senses, master_lb, TN, P)
    Gamma <- master_result$Gamma
    quadratic_penalty <- max(master_result$quadratic_penalty, 0)
    linear_penalty <- max(master_result$linear_penalty, 0)
    for(t in 1:TN)
    {
      eta[[t]] <- regularize_eta(master_result$eta[[t]], parameters$epsilon)
    }
    
    kappa <- vector("list", TN)
    
    for (t in 1:TN) {
      Keta <- calculate_Keta(Km[[t]], eta[[t]])
      svm_result <- solve_classification_svm(Keta, y[[t]], parameters$C, parameters$epsilon)
      
      alpha[[t]] <- svm_result$alpha #regularized alpha (may be infeasible)
      b[[t]] <- svm_result$b
      
      alpha_cut <- svm_result$alpha_original
      J[t] <- svm_result$objective_original
      
      #model$alpha_original is already multiplied by y
      alpha_sum <- t(alpha_cut) %*% y[[t]]
      kappa[[t]] <- calucalte_kappa(alpha_cut, Km[[t]], P)
      
      #add constraint Gamma_t + Eta_t*kappa >= alpha_sum
      constraint_alpha <- numeric(NVars)
      constraint_alpha[t] <- 1
      constraint_alpha[(TN+1+(t-1)*P):(TN+1+(t)*P-1)] <- kappa[[t]]
      
      master_constraint_matrix <- rbind(master_constraint_matrix, constraint_alpha)
      master_rhs <- c(master_rhs, alpha_sum)
      master_senses <- c(master_senses, "G")  
      
      new_svm_constraints <- rbind(new_svm_constraints, constraint_alpha)
      new_svm_rhs <- c(new_svm_rhs, alpha_sum)
    }
    
    obj <- sum(J) + linear_penalty + quadratic_penalty
    if(obj < UB)
    {
      UB <- obj
      incumbent_solution = list(eta = eta, alpha = alpha, b = b)
    }
    
    LB <- sum(Gamma) + linear_penalty + quadratic_penalty
    
    Gap <- (UB-LB)/max(abs(LB),parameters$epsilon)
    
    objectives <- c(objectives, obj)
    
    info <- c(obj, UB, LB, Gap, sum(Gamma), sum(J), quadratic_penalty, linear_penalty)
    iteration_info <- rbind(iteration_info, info)
    
    if(iteration %% 5 == 1)
      print(sprintf("%d %f %f %f %f%% %f %f %f %f",iteration, obj, UB, LB, Gap*100, sum(Gamma), sum(J), quadratic_penalty, linear_penalty))
    
    if(Gap <= parameters$optimality_gap | iteration >= parameters$iteration_count)
    {
      terminate <- TRUE
    }
    
    if(terminate & iteration %% 5 != 1)
      print(sprintf("%d %f %f %f %f%% %f %f %f %f",iteration, obj, UB, LB, Gap*100, sum(Gamma), sum(J), quadratic_penalty, linear_penalty))
    
    iteration <- iteration + 1
  }
  
  state <- list(alpha = incumbent_solution$alpha, b = incumbent_solution$b, eta = incumbent_solution$eta, 
                objectives = objectives, iteration_info = iteration_info, parameters = parameters)
  
  new_cuts_info = list(SVM_cuts = list(matrix = new_svm_constraints, rhs = new_svm_rhs))
  
  output <- list(state = state, new_cuts_info = new_cuts_info)
  
  return(output)
}
