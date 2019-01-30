library(Rcplex)

solve_classification_svm <- function(K, y, C, epsilon) {
  N <- length(y)
  yyK <- (y %*% t(y)) * K
  
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  result <- Rcplex(cvec = rep(-1, N), 
                   Amat = matrix(y, nrow = 1, byrow = TRUE), 
                   bvec = 0, 
                   Qmat = yyK, 
                   lb = rep(0, N),
                   ub = rep(C, N),
                   control = opts,
                   objsense = "min",
                   sense = "E")
  
  alpha <- result$xopt[1:N]
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
  
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  Qmat <- NULL
  if(!is.null(obj_quadratic))
    Qmat <- 2*obj_quadratic
  
  result <- Rcplex(cvec = obj_coef, 
                   Qmat = Qmat,
                   Amat = constraints_matrix, 
                   bvec = rhs, 
                   lb = lb,
                   control = opts,
                   objsense = "min",
                   sense = senses)
  
  #order of variables: Gamma_1, ..., Gamma_T, eta_1_1,...,eta_1_P,eta_2_1,...,eta_2_P,...,eta_T_1,...,eta_T_P
  sln <- result$xopt
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
    quadratic_penalty <- (t(sln) %*% obj_quadratic) %*% (sln)
  
  linear_penalty <- sum(sln*obj_coef) - sum(Gamma)
  
  objective <- result$obj
  
  output <- list(Gamma = Gamma, quadratic_penalty = quadratic_penalty, linear_penalty = linear_penalty, eta = eta, objective = objective)
  return(output)
}