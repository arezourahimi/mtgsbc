pdist <- function(X1, X2) {
  if (identical(X1, X2) == TRUE) {
    D <- as.matrix(dist(X1))
  }
  else {
    D <- as.matrix(dist(rbind(X1, X2)))
    D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
  }
  return(D)
}

read_pathways <- function(name) {
  symbols_lines <- read.table(sprintf("msigdb/%s.gmt", name), header = FALSE, sep = ",", stringsAsFactor = FALSE)
  pathways <- vector("list", nrow(symbols_lines))
  for (line in 1:nrow(symbols_lines)) {
    symbols_entries <- strsplit(symbols_lines[line, 1], "\t")
    pathways[[line]]$name <- symbols_entries[[1]][1]
    pathways[[line]]$link <- symbols_entries[[1]][2]
    pathways[[line]]$symbols <- sort(symbols_entries[[1]][-2:-1])
  }
  return(pathways)
}

calculate_Keta <- function(Km, eta) {
  P <- dim(Km)[3]
  
  Keta <- eta[1] * Km[,,1]
  for (m in 2:P) {
    Keta <- Keta + eta[m] * Km[,,m]
  }
  
  return(Keta)
}

#Calculates the penalty info for a given list of eta vectors, and a corresponding penalty matrix
#Returns: penalty value, gradient, and inner product of gradiennt and eta
calculate_penalty_info <- function(eta, penalty_matrix)
{
  TN <- length(eta)
  NV <- nrow(penalty_matrix)
  
  eta_linear = c()
  for(s in 1:TN)
  {
    eta_linear <- c(eta_linear, eta[[s]])
  }
  
  epm = t(eta_linear) %*% penalty_matrix
  
  penalty_info = list()
  penalty_info$gradient = 2*epm
  penalty_info$penalty = epm %*% eta_linear
  penalty_info$gradient_product = 2*penalty_info$penalty
  
  return(penalty_info)
}

#Calculates the coefficients of kernel weights in the objective function of dual SVM problem for a given alpha
calucalte_kappa <- function(alpha, Km, P)
{
  kappa <- numeric(P)
  
  for(m in 1:P)
  {
    kappa[m] <- 0.5*t(alpha) %*% Km[,,m] %*% alpha
  }
  return(kappa)
}

#Removes the numerically insignificat kernel weights
regularize_eta <- function(eta, epsilon)
{
  et <- eta / sum(eta)
  et[et < epsilon] <- 0
  et <- et / sum(et)
  return(et)
}

#Given a matrix of grouping weights of T types (either corresponding to coclusters, or any other possible values),
#populates a (TN*P) by (TN*P) matrix equivelent to the quadratic penalty term. 
generate_penalty_matrix <- function(grouping_weight_matrix, P)
{
  TN <- nrow(grouping_weight_matrix)
  penalty_matrix <- matrix(0, nrow = P*TN, ncol = P*TN)
  grouping_weight <- colSums(grouping_weight_matrix)-diag(grouping_weight_matrix)
  
  
  #Note: the (sm,qn) entry of penalty_matrix is equal to zero if n!=m
  for(m in 1:P)
  {
    for(s in 1:TN)
    {
      sm = (s-1)*P+m
      
      for(q in 1:TN)
      {
        if(s==q)
        {
          penalty_matrix[sm,sm] <- 2*grouping_weight[s]
        }else
        {
          qm = (q-1)*P+m
          penalty_matrix[sm,qm] <- -2*grouping_weight_matrix[s,q]
        }
      }
    }
  }
  
  return(penalty_matrix)
}

#Converts a square matrix to its sparse linear counterpart, used by MOSEK.
linear_index_penalty_matrix <- function(penalty_matrix, shift = 0)
{
  N <- nrow(penalty_matrix)
  I <- matrix(1:N, N, N, byrow = FALSE)
  J <- matrix(1:N, N, N, byrow = TRUE)
  lt <- lower.tri(I, diag = FALSE)
  
  i1 <- I[lt]
  j1 <- J[lt]
  v1 <- penalty_matrix[lt]
  
  i2 <- 1:N
  j2 <- 1:N
  v2 <- diag(penalty_matrix)
  
  i <- c(i1,i2)
  j <- c(j1,j2)
  v <- c(v1, v2)
  
  nonzero <- which(v!=0)
  
  q <- list(i=i[nonzero]+shift, j=j[nonzero]+shift, v = v[nonzero])
  
  return(q)
}

#Given eta (in matrix form, with T columns and P rows), finds a pool of diverse near optimal solutions using a variation of GA
generate_block_clusters <- function(eta, K, lambda_ratio = 0, min_cluster_size = 2,  diversity_dist = 10, required_diverse_solutions = 10)
{
  diverse_solutions <- coclustering_heuristic_algorithm(eta = eta, K = K, lambda_1 = 1, lambda_2 = lambda_ratio,
                                                         min_group_size = min_cluster_size, max_group_size_ratio = 1,
                                                         mutations = 2, swaps = 1,
                                                         max_local_nonimproving = 10000,
                                                         max_global_nonimproving = 5, 
                                                         diversity_dist = diversity_dist,
                                                         required_diverse_solutions = required_diverse_solutions,
                                                         reoptimize = TRUE)
  cc <- diverse_solutions[[1]] #global best solution
  return(list(cohorts = cc[[1]], pathways = cc[[2]], diverse_solutions = diverse_solutions))  
}

generate_block_penalty_info_from_coclusters <- function(block_clusters)
{
  grouping_matrices <- list()
  grouping_matrices$cohorts <- calculate_cocluster_matrix(block_clusters$cohorts)
  grouping_matrices$pathway_cohorts <- calculate_cocluster_matrix(block_clusters$pathways, block_clusters$cohorts)
  
  return(generate_block_penalty_info(grouping_matrices))
}

generate_block_penalty_info <- function(grouping_matrices)
{
  #grouping_matrices contains these two matrices: cohorts & pathway_cohorts
  
  TN <- ncol(grouping_matrices$pathway_cohorts)
  P <- nrow(grouping_matrices$pathway_cohorts)
  
  penalty_matrix <- generate_penalty_matrix(grouping_matrices$cohorts, P)
  
  penalty_vector <- 1-as.vector(grouping_matrices$pathway_cohorts)
  
  return(list(penalty_matrix = penalty_matrix, penalty_vector = penalty_vector))
}

#Given two group_index vectors, generates a matrix with element (i,j) equal to 1 
#if the i'th element of row elements and the j'th element of column elements belong to the same cluster
#if col_groups is NULL, then the two index vectors are assumed to be identical, which results in a square matrix
calculate_cocluster_matrix <- function(row_group_indexes, col_group_indexes = NULL)
{
  if(is.null(col_group_indexes))
    col_group_indexes <- row_group_indexes
  
  N1 <- length(row_group_indexes)
  N2 <- length(col_group_indexes)
  cocluster_matrix <- matrix(0, nrow = N1, ncol = N2)
  for(i in 1:N1)
  {
    for(j in 1:N2)
    {
      if(row_group_indexes[i] == col_group_indexes[j])
      {
        cocluster_matrix[i,j] <- 1
      }
    }
  }
  
  return(cocluster_matrix)
}

#Calculates the disatnce of two clusters, defined as the minimum number of displaced elements over K! permutations of cluster labels
calculate_cluster_distance <- function(g1, g2, K)
{
  clusters1 <- list()
  clusters2 <- list()
  
  for(k in 1:K)
  {
    clusters1[[k]] <- which(g1 == k)
    clusters2[[k]] <- which(g2 == k)
  }
  
  perms <- permutations(n=K)
  
  min_dist = Inf
  for(p in 1:nrow(perms))
  {
    cls2 <- clusters2[perms[p,]]
    dist <- 0
    for(k in 1:K)
    {
      intersection <- intersect(clusters1[[k]], cls2[[k]])
      union <- union(clusters1[[k]], cls2[[k]])
      dist <- dist + length(setdiff(union, intersection))
    }
    
    if(dist < min_dist)
      min_dist <- dist
  }
  
  return(min_dist)
}

permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else 
  {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

#Yields different combinations of C, Lambda_1, and Ratio (gamma) for the coclustering step
get_cross_validation_tuples_coclustering_step <- function()
{
  C_set <- 10^(-4:5)
  lambda_set <- 10^(-1:4)
  lambda_ratio_range <- 5^(-1:2)
  
  tuples <- as.matrix(expand.grid(C_set, lambda_set, lambda_ratio_range))
  colnames(tuples) <- c("C", "Lambda", "Ratio")
  
  return(tuples)
}

#For the final MTMKL step, we use the same tuples used for the coclustering step
get_cross_validation_tuples_mtmkl <- function()
{
  C_set <- 10^(-4:5)
  lambda_set <- 10^(-1:4)
  lambda_ratio_range <- 5^(-1:2)
  
  tuples <- as.matrix(expand.grid(C_set, lambda_set, lambda_ratio_range))
  colnames(tuples) <- c("C", "Lambda", "Ratio")
  
  return(tuples)
}
