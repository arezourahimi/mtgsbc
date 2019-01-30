library(Rmosek)

solve_classification_svm <- function(K, y, C, epsilon) {
  N <- length(y)
  yyK <- (y %*% t(y)) * K
  
  problem <- list()
  problem$sense <- "min"
  problem$c <- rep(-1, N)
  problem$A <- Matrix(y, nrow = 1, byrow = TRUE, sparse = TRUE)
  problem$bc <- rbind(blc = 0, buc = 0) 
  
  bux <- rep(C, N)
  
  problem$bx <- rbind(blx = rep(0, N), bux = bux)
  
  I <- matrix(1:N, N, N, byrow = FALSE)
  J <- matrix(1:N, N, N, byrow = TRUE)
  problem$qobj <- list(i = I[lower.tri(I, diag = TRUE)],
                       j = J[lower.tri(J, diag = TRUE)],
                       v = yyK[lower.tri(yyK, diag = TRUE)])
  
  problem$iparam <- list(INTPNT_MULTI_THREAD=0, NUM_THREADS=1)
  
  
  opts <- list()
  opts$verbose <- 0
  result <- mosek(problem, opts)
  
  alpha <- result$sol$itr$xx[1:N]
  alpha_original <- alpha
  alpha[alpha < +C * epsilon] <- 0
  alpha[alpha > +C * (1 - epsilon)] <- +C
  
  objective <- sum(alpha) - 0.5 * (t(alpha) %*% yyK) %*% (alpha)
  objective <- objective * (objective >= 0)
  
  objective_original <- sum(alpha_original) - 0.5 * (t(alpha_original) %*% yyK) %*% (alpha_original)
  objective_original <- objective_original * (objective_original >= 0)
  
  support_indices <- which(alpha != 0)
  active_indices <- which(alpha != 0 & alpha < C)
  
  if (length(active_indices) > 0) {
    b <- mean(y[active_indices] * (1 - yyK[active_indices, support_indices] %*% alpha[support_indices]))
  } else {
    b <- 0
  }
  
  model <- list(alpha = alpha * y, b = b, objective = objective, alpha_original = alpha_original*y, objective_original = objective_original)
  
  return(model)
}

solve_multitask_quadratic_master_problem <- function(obj_coef, obj_quadratic, constraints_matrix, 
                                                     rhs, senses, lb, TN, P) {
  blc <- rhs
  buc <- rhs
  blc[which(senses=="L")] <- -Inf
  buc[which(senses=="G")] <- Inf
  
  problem <- list()
  problem$sense <- "min"
  problem$c <- obj_coef
  problem$A <- Matrix(constraints_matrix)
  problem$bc <- rbind(blc = blc, buc = buc) 
  problem$bx <- rbind(blx = lb, bux = Inf)
  
  problem$iparam <- list(INTPNT_MULTI_THREAD=0, NUM_THREADS=1)
  
  if(!is.null(obj_quadratic))
  {
    if(class(obj_quadratic) == "matrix")
    {
      #obj_quadratic already takes the initial TN+1 variables so no need for shifting
      #obj_quadratic is in the original form of values. In the follwoing function
      #main diagonal elements are multiplied by 2
      problem$qobj <- linear_index_penalty_matrix(obj_quadratic)
    }else #obj_quadratic is already in the linear format
    {
      problem$qobj <- obj_quadratic
    }
  }
  
  opts <- list()
  opts$verbose <- 0
  result <- mosek(problem, opts)
  
  #order of variables: Gamma_1, ..., Gamma_T, rho, eta_1_1,...,eta_1_P,eta_2_1,...,eta_2_P,...,eta_T_1,...,eta_T_P
  sln <- result$sol$itr$xx
  Gamma <- sln[1:TN]
  
  i <- TN+1
  eta <- vector("list", TN)
  for(t in 1:TN)
  {
    eta[[t]] <- sln[i:(i+P-1)]
    i <- i+P
  }
  
  quadratic_penalty <- 0
  if(!is.null(obj_quadratic))
  {
    if(class(obj_quadratic) == "matrix")
    {
      quadratic_penalty <- (t(sln) %*% obj_quadratic) %*% (sln)
    }else
    {
      for (li in 1:length(obj_quadratic$v)) {
        ms <- obj_quadratic$i[li]
        mq <- obj_quadratic$j[li]
        if(ms == mq)
        {
          v <- 0.5*obj_quadratic$v[li]
        }else
        {
          v <- obj_quadratic$v[li]
        }
        quadratic_penalty <- quadratic_penalty + sln[ms]*sln[mq]*v
      }
    }
  }
  
  linear_penalty <- sum(sln*obj_coef) - sum(Gamma)
  objective <- sum(Gamma) + linear_penalty + quadratic_penalty
  
  output <- list(Gamma = Gamma, quadratic_penalty = quadratic_penalty, linear_penalty = linear_penalty, eta = eta, objective = objective)
  return(output)
}