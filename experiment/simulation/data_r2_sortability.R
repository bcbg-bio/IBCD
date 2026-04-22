library(inspre)
library(tidyverse)
library(hdf5r)
library(progress)
library(mc2d)
library(doParallel)
library(foreach)
library(progressr)
library(doRNG)
library(tidyverse)
library(mc2d)

multiple_iv_reg_UV <- function(target, .X, .targets){
  inst_obs <- .targets == target
  control_obs <- .targets == "control"
  n_target <- sum(inst_obs)
  n_control <- sum(control_obs)
  n_total <- n_target + n_control
  X_target <- .X[inst_obs, ]
  X_control <- .X[control_obs, ]
  this_feature <- which(colnames(.X) == target)
  X_exp_target <- X_target[, this_feature]
  X_exp_control <- X_control[, this_feature]
  beta_inst_obs <- sum(X_exp_target)/n_target
  sums_target <- colSums(X_target)
  beta_hat <- sums_target/sums_target[this_feature]
  resid <- rbind(X_target - outer(X_exp_target, beta_hat), X_control - outer(X_exp_control, beta_hat))
  V_hat <- (t(resid)%*%resid)/nrow(resid)
  se_hat <- sqrt(diag(V_hat)/(n_target*beta_inst_obs**2))
  names(beta_hat) <- paste(names(beta_hat), "beta_hat", sep="_")
  names(se_hat) <- paste(names(se_hat), "se_hat", sep="_")
  return(list(
    beta_se = as.data.frame(c(list(target = target, beta_obs = beta_inst_obs), 
                              as.list(c(beta_hat, se_hat)))),
    U_i=1/(n_target*beta_inst_obs**2),
    V_hat=V_hat
  ))
  # return(beta_se=as.data.frame(c(list(target = target, beta_obs = beta_inst_obs), as.list(c(beta_hat, se_hat)))), U_i=1/n_target*beta_inst_obs**2, V_hat=V_hat)
}

calc_R2_i <- function(covY, i){
  covY_mi <- covY[-i, -i]
  covY_i <- covY[i, -i]
  beta_hat <- solve(covY_mi, covY_i)
  R2 <- sum(beta_hat * colSums(covY_mi * beta_hat))
  return(R2)
}


calc_R2 <- function(G){
  D <- nrow(G)
  H <- solve(diag(D) - G)
  covY <- t(H) %*% H
  
  R2s <- purrr::map_dbl(1:nrow(G), ~ calc_R2_i(covY, .x))
  return(R2s)
}


build_design <- function(M){
  D <- nrow(M)
  N <- D*(D+1)/2
  B <-matrix(0, D, N)
  
  k = 0
  for(i in 1:D){
    l = 0
    for(j in i:D){
      k = k+1
      B[i, k] = M[i, j]
      if(i!=j){
        B[i+l, k] = M[i, j]
      }
      l = l+1
    }
  }
  return(B)
}


generate_data_inhibition <- function(G, N_cont, N_int, int_beta=-2, noise='gaussian', normalize=TRUE, net_vars=NULL, min_v=NULL){
  D <- nrow(G)
  int_sizes <- rep(N_int, D) # rnbinom(D, mu=N_int - N_int/10, size=N_int/10) + N_int/10
  
  N_int_cs <- cumsum(c(0, int_sizes))
  XB <- matrix(0, nrow=N_int_cs[length(N_int_cs)], ncol=D)
  
  #for(d in 1:D){
  #  start = N_int_cs[d]
  #  end = N_int_cs[d+1]
  #  XB[(start+1):end, d] = 1
  #}
  
  # observation only
  if (N_int > 0) {
    for(d in 1:D){
      start = N_int_cs[d]
      end = N_int_cs[d+1]
      if (start < end) {
        XB[(start+1):end, d] = 1
      }
    }
  }
  
  if(is.null(net_vars)){
    net_vars <- colSums(G**2)
    eps_vars <- max(0.9, max(net_vars)) - net_vars + 0.1
  } else{
    net_vars_obs <- colSums(G**2)
    rescale <- sqrt(net_vars/net_vars_obs)
    rescale[is.infinite(rescale)] = 1
    G <- t((t(G)*rescale))
    
    G_eff <- colSums(G**2) + rowSums(G**2)
    A <- sqrt(outer(1/G_eff, 1/G_eff))
    GA_int <- G*A
    GA_int_eff <- colSums(GA_int**2) + rowSums(GA_int**2)
    
    print(summary(abs(G)[abs(G) > 0.0001]))
    G <- GA_int * sqrt(mean(G_eff)/mean(GA_int_eff))
    if(!is.null(min_v)){ G[abs(G) < min_v] <- 0 }
    print(summary(abs(G)[abs(G) > 0.0001]))
    
    net_vars <- colSums(G**2)
    eps_vars <- max(0.9, max(net_vars)) - net_vars + 0.1
  }
  total_var <- net_vars + eps_vars
  
  if(noise == 'gaussian'){
    eps_cont <- t(matrix(rnorm(D*N_cont, sd=sqrt(eps_vars)), nrow=D, ncol=N_cont))
    eps_int <- t(matrix(rnorm(D*sum(int_sizes), sd=sqrt(eps_vars)), nrow=D, ncol=sum(int_sizes)))
  } else{
    stop('NotImplementedError')
  }
  Y_cont <- eps_cont %*% solve(diag(D) - G)
  if(normalize){
    # This mimics perturb-seq normalization
    mu_cont <- colMeans(Y_cont)
    sd_cont <- apply(Y_cont, 2, sd)
    
    Y_cont <- t((t(Y_cont) - mu_cont)/sd_cont)
    Y_int <- (t(t(XB) * int_beta * sd_cont) + eps_int) %*% solve(diag(D) - G)
    Y_int <- t((t(Y_int) - mu_cont)/sd_cont)
    
    # Also need to normalize the graph
    R <- get_tce(get_observed(G), normalize=sd_cont)
    G <- get_direct(R)$G
    print(summary(abs(G)[abs(G) > 0.0001]))
    # TODO: update eps_vars
  } else{
    Y_int <- (t(t(XB) * int_beta) + eps_int) %*% solve(diag(D) - G)
    R <- get_tce(get_observed(G))
  }
  Y <- rbind(Y_cont, Y_int)
  
  colnames(Y) <- paste0("V", 1:D)
  targets <- c(rep("control", N_cont), paste0("V", rep(1:D, times=int_sizes)))
  return(list(Y = Y, targets = targets, G = G, R = R, int_beta=int_beta, eps_vars=eps_vars))
}

generate_dataset <- function(D, N_cont, N_int, int_beta=-2,
                             graph = 'scalefree', v = 0.2, p = 0.4,
                             DAG = FALSE, C = floor(0.1*D), noise = 'gaussian',
                             model = 'inhibition', net_vars = NULL, min_v=NULL){
  if(!is.null(net_vars)){
    net_vars = rep(net_vars, D)
  }
  G <- generate_network(D, graph, p, v, DAG)
  if(C > 0){
    if(graph == "scalefree"){
      new_vars <- matrix(mc2d::rpert(D*C, 0.01, v^(log(D)/log(log(D))), 2*v), nrow=C)
    } else if(graph == "random"){
      new_vars <- matrix(mc2d::rpert(D*C, 0.01, v^(log(D)), 2*v), nrow=C)
    }
    
    G <- cbind(rbind(G, new_vars), matrix(0, nrow=D+C, ncol=C))
  }
  if(model == 'inhibition'){
    data = generate_data_inhibition(G, N_cont, N_int, int_beta, noise, net_vars=net_vars, min_v=min_v)
  } else if (model == 'knockout'){
    data = generate_data_knockout(G, N_cont, N_int, noise)
  }
  
  Y <- data$Y
  G <- data$G
  var_all <- apply(Y %*% G, 2, var)[1:D]
  var_obs <- apply(Y[,1:D] %*% G[1:D, 1:D], 2, var)
  if(C > 0){
    var_conf <- apply(Y[,(D+1):(D+C)] %*% G[(D+1):(D+C), 1:D], 2, var)
  } else {
    var_conf <- rep(0, length=D)
  }
  var_eps <- 1-var_all
  
  return(list(Y=Y[1:(N_cont+N_int*D), 1:D], targets=data$targets[1:(N_cont+N_int*D)], R=data$R[1:D, 1:D],
              G=G[1:D, 1:D], var_all=var_all, var_obs=var_obs,
              var_conf=var_conf, var_eps=var_eps, int_beta=data$int_beta, eps_vars=data$eps_vars))
}


int_beta_default <- -2 
v_default <- 0.25 
graph_type <- "scalefree" 
DAG_flag <- TRUE
C_default <- 0
min_v <- v_default/2
replicate_seeds <- 42:51 

run_one <- function(D, N_int, p, seed, out_dir,
                    int_beta = int_beta_default,
                    v = v_default,
                    graph = graph_type,
                    DAG = DAG_flag,
                    C = C_default,
                    net_vars = net_vars,
                    min_v = min_v) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("[D=%d | N_int=%d | p=%.3f | seed=%d] -> %s",
                  D, N_int, p, seed, out_dir))
  set.seed(seed)
  
  N_cont <- D * N_int
  
  # --- Generate dataset ---
  dataset <- generate_dataset(
    D = D,
    N_cont = N_cont,
    N_int = N_int,
    int_beta = int_beta,
    graph = graph,
    v = v,
    p = p,
    DAG = DAG,
    C = C,
    net_vars = net_vars,
    min_v = min_v
  )
  
  G <- dataset$G
  R <- dataset$R
  Y_fixed <- dataset$Y
  targets_fixed <- dataset$targets
  writeLines(targets_fixed, file.path(out_dir, "targets.txt"))
  
  Y_C <- Y_fixed[targets_fixed == "control", ]
  Sigma <- t(Y_C) %*% Y_C / nrow(Y_C)
  xi <- Sigma^2
  diag(xi) <- 0
  xi_norm <- xi / sqrt(sum(xi^2))
  
  write.csv(G, file = file.path(out_dir, "G_matrix.csv"), row.names = FALSE)
  write.csv(R, file = file.path(out_dir, "R_matrix.csv"), row.names = FALSE)
  write.csv(Y_fixed, file = file.path(out_dir, "Y_matrix.csv"), row.names = FALSE)
  write.csv(xi_norm, file = file.path(out_dir, "xi_true.csv"), row.names = FALSE)
  
  genes <- paste0("V", 1L:D)
  U_diag <- numeric(D)
  V_sum <- matrix(0, D, D)
  R_hat <- matrix(0, D, D)
  SE_hat <- matrix(0, D, D)
  
  for (i in seq_along(genes)) {
    gene <- genes[i]
    message("Running IV for target: ", gene)
    res <- multiple_iv_reg_UV(gene, Y_fixed, targets_fixed)
    U_diag[i] <- res$U_i
    V_sum <- V_sum + res$V_hat
    R_hat[i, ] <- unlist(res$beta_se[1, grep("_beta_hat$", names(res$beta_se))])
    SE_hat[i, ] <- unlist(res$beta_se[1, grep("_se_hat$", names(res$beta_se))])
  }
  
  U <- diag(U_diag)
  V <- V_sum / D
  write.csv(R_hat, file = file.path(out_dir, "R_true.csv"), row.names = FALSE)
  write.csv(SE_hat, file = file.path(out_dir, "SE_hat.csv"), row.names = FALSE)
  write.csv(U,      file = file.path(out_dir, "U_true.csv"),    row.names = FALSE)
  write.csv(V,      file = file.path(out_dir, "V_true.csv"),    row.names = FALSE)
}

#######################################
############ int only  ################
#######################################

base_dir <- file.path(getwd(), "int_r2")

if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}

D_50 <- 50
p_50 <- 0.066
subfolders_50d <- c(5, 15, 25, 50, 75, 100)
net_vars <- 0.4

for (Nint in subfolders_50d) {
  for (seed in replicate_seeds) {
    out_dir <- file.path(base_dir, "50d", as.character(Nint), "sf", as.character(seed))
    run_one(D = D_50, N_int = Nint, p = p_50, seed = seed, out_dir = out_dir, net_vars = net_vars, min_v = min_v)
  }
}

# -----------------------------
D_250 <- 250
Nint_250 <- 100
p_250 <- 0.124
net_vars <- 0.5

for (seed in replicate_seeds) {
  out_dir <- file.path(base_dir, "250d", "sf", as.character(seed))
  run_one(D = D_250, N_int = Nint_250, p = p_250, seed = seed, out_dir = out_dir, net_vars = net_vars, min_v = min_v)
}



#######################################
############ obs only  ################
#######################################


int_beta_default <- -2
v_default <- 0.25
graph_type <- "scalefree" 
DAG_flag <- TRUE
C_default <- 0
N_int_default <- 0
min_v <- v_default/2
replicate_seeds <- 42:51  # 10 replicates

run_one <- function(D, variable, N_int = N_int_default, p, seed, out_dir,
                    int_beta = int_beta_default,
                    v = v_default,
                    graph = graph_type,
                    DAG = DAG_flag,
                    C = C_default,
                    net_vars = net_vars,
                    min_v = min_v) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("[D=%d | N_int=%d | p=%.3f | seed=%d] -> %s",
                  D, N_int, p, seed, out_dir))
  set.seed(seed)
  
  N_cont <- 2 * D * variable
  
  # --- Generate dataset ---
  dataset <- generate_dataset(
    D = D,
    N_cont = N_cont,
    N_int = N_int,
    int_beta = int_beta,
    graph = graph,
    v = v,
    p = p,
    DAG = DAG,
    C = C,
    net_vars = net_vars,
    min_v = min_v
  )
  
  G <- dataset$G
  R <- dataset$R
  Y_fixed <- dataset$Y
  targets_fixed <- dataset$targets
  writeLines(targets_fixed, file.path(out_dir, "targets.txt"))
  
  Y_C <- Y_fixed[targets_fixed == "control", ]
  Sigma <- t(Y_C) %*% Y_C / nrow(Y_C)
  xi <- Sigma^2
  diag(xi) <- 0
  xi_norm <- xi / sqrt(sum(xi^2))
  
  write.csv(G, file = file.path(out_dir, "G_matrix.csv"), row.names = FALSE)
  write.csv(R, file = file.path(out_dir, "R_matrix.csv"), row.names = FALSE)
  write.csv(Y_fixed, file = file.path(out_dir, "Y_matrix.csv"), row.names = FALSE)
  write.csv(xi_norm, file = file.path(out_dir, "xi_true.csv"), row.names = FALSE)
  
}

base_dir <- file.path(getwd(), "obs_r2")

if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}


D_50 <- 50
p_50 <- 0.066
variables_50d <- c(5, 15, 25, 50, 75, 100)
net_vars <- 0.4

for (variable in variables_50d) {
  for (seed in replicate_seeds) {
    out_dir <- file.path(base_dir, "50d", as.character(variable), "sf", as.character(seed))
    run_one(D = D_50, variable = variable, p = p_50, seed = seed, out_dir = out_dir, net_vars = net_vars, min_v = min_v)
  }
}

# -----------------------------
D_250 <- 250
p_250 <- 0.124
variable <- 100
net_vars <- 0.5

for (seed in replicate_seeds) {
  out_dir <- file.path(base_dir, "250d", "sf", as.character(seed))
  run_one(D = D_250, variable = variable, p = p_250, seed = seed, out_dir = out_dir, net_vars = net_vars, min_v = min_v)
}



