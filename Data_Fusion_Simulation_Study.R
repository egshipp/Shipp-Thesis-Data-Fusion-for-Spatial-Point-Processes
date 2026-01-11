# Required Packages -------------------------------------------------------------
library(spatstat)
library(fields)
library(FastGP)
library(MASS)
library(ggplot2)
library(coda)

# True LGCP --------------------------------------------------------------------
## Simulate Covariate X(s)
win <- owin(xrange = c(0, 10), yrange = c(0, 10))

grid_res <- 10

cell_size <- diff(win$xrange) / grid_res

x_seq <- seq(win$xrange[1] + cell_size/2,
             win$xrange[2] - cell_size/2,
             by = cell_size)
y_seq <- seq(win$yrange[1] + cell_size/2,
             win$yrange[2] - cell_size/2,
             by = cell_size)

coords <- as.matrix(expand.grid(y_seq, x_seq))
grid_coords <- expand.grid(x = y_seq, y = x_seq)

dists <- as.matrix(dist(coords))
n <- nrow(coords)

S <- 0.2 * exp(-dists/4)

mu <- rep(0, n)

covariate <- as.vector(rcpp_rmvnorm(1,S,mu))
# save(covariate, file = "sim_covariate.Rdata")

cov_field <- matrix(covariate - mean(covariate),
                    nrow = grid_res,
                    ncol = grid_res,
                    byrow = TRUE)

cov_field_ppp <- ppp(x = grid_coords$x,
                     y = grid_coords$y,
                     window = win,
                     marks = covariate)

## Simulate Gaussian random field z(s)
sigma_2 <- 0.10

S_z <-  sigma_2 * exp(-dists/4)

z <- as.vector(rcpp_rmvnorm(1,S_z,mu))
# save(z, file = "z.RData")

z_mat <- matrix(z, nrow = length(x_seq), ncol = length(y_seq))

z_ppp <- ppp(x = grid_coords$x,
             y = grid_coords$y,
             window = win,
             marks = z)

## Simulate entire true LGCP
b_0 <- 1.5
b_1 <- 3

lambda <- exp(b_0 + b_1*(cov_field) + z)
# save(lambda, file = "sim_lambda.RData")
lambda_im <- im(lambda, xcol = x_seq, yrow = y_seq)

lgcp_sim <- rpoispp(lambda_im)
# save(lgcp_sim, file = "lgcp_sim.RData")

## Discretize true point process using spatstat
lgcp_discretize <- pixellate(lgcp_sim, eps = 1)

# Source 1  --------------------------------------------------------------------

## Creating sub-region for Source 1 domain
nrow <- length(lgcp_discretize$yrow)
ncol <- length(lgcp_discretize$xcol)

x_min_subwindow1 <- 0
x_max_subwindow1 <- 5
y_min_subwindow1 <- 5
y_max_subwindow1 <- 10

x_min_subwindow2 <- 5
x_max_subwindow2 <- 10
y_min_subwindow2 <- 0
y_max_subwindow2 <- 5

sub_window1 <- owin(xrange = c(x_min_subwindow1, x_max_subwindow1), yrange = c(y_min_subwindow1, y_max_subwindow1))

sub_window2 <- owin(xrange = c(x_min_subwindow2, x_max_subwindow2), yrange = c(y_min_subwindow2, y_max_subwindow2))

## Simulating Source 1 point pattern
lambda_1 <- lgcp_discretize

lgcp_1_sub1 <- rpoispp(lambda_1[sub_window1])
lgcp_1_sub2 <- rpoispp(lambda_1[sub_window2])

lgcp_1 <- superimpose(lgcp_1_sub1, lgcp_1_sub2, W = owin(c(0,10), c(0,10)))
# save(lgcp_1, file = "lgcp_1.RData")

# Source 2 ----------------------------------------------------------------------
tau_2 <- 0.4
S_g <- tau_2 * exp(-dists/4)
alpha <- -0.2

g <- as.vector(rcpp_rmvnorm(1,S_g,alpha))
# save(g, file = "g.RData")

exp_g <- exp(g)

lambda_2 <- lgcp_discretize * exp_g

lgcp_2 <- rpoispp(lambda_2)
# save(lgcp_2, file = "lgcp_2.RData")

# Data frame creation ----------------------------------------------------------

## X grid
X_grid <- as.data.frame(coords)
colnames(X_grid) <- c("x", "y")
X_grid$covariate <- covariate

## Source 1 (small var)
X_1 <- as.data.frame(lgcp_1)
nn_X_1 <- nncross(lgcp_1, cov_field_ppp)
X_1$covariate <- cov_field_ppp$marks[nn_X_1$which]

# making a mask for X-grid to be used with source 1 in log-likelihood function
inside_sub1 <- with(X_grid,
                    x >= x_min_subwindow1 & x <= x_max_subwindow1 &
                      y >= y_min_subwindow1 & y <= y_max_subwindow1)

inside_sub2 <- with(X_grid,
                    x >= x_min_subwindow2 & x <= x_max_subwindow2 &
                      y >= y_min_subwindow2 & y <= y_max_subwindow2)

X_grid$mask_source1 <- inside_sub1 | inside_sub2

## Source 2 (large var)
X_2 <- as.data.frame(lgcp_2)
nn_X_2 <- nncross(lgcp_2, cov_field_ppp)
X_2$covariate <- cov_field_ppp$marks[nn_X_2$which]

# MCMC --------------------------------------------------------------------------

## log likelihood function
loglike <- function(parameters, data) {
  
  idx_mask1 <- which(data$X_grid$mask_source1)
  
  log_lambda_points1 <- parameters$beta[1] + parameters$beta[2] * data$X_1$covariate + parameters$z[data$nn_index_1]
  term1 <- sum(log_lambda_points1)
  
  log_lambda_grid1_full <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + parameters$z
  
  lambda_grid1_masked <- exp(log_lambda_grid1_full[idx_mask1])
  term2 <- sum(lambda_grid1_masked * data$cell_area)
  
  log_lambda_points2 <- parameters$beta[1] + parameters$beta[2] * data$X_2$covariate + parameters$g[data$nn_index_2] + parameters$z[data$nn_index_2]
  term3 <- sum(log_lambda_points2)
  
  log_lambda_grid2 <- parameters$beta[1] + parameters$beta[2] * data$X_grid$covariate + parameters$g + parameters$z
  lambda_grid2 <- exp(log_lambda_grid2)
  term4 <- sum(lambda_grid2 * data$cell_area)
  
  likelihood <- (term1 - term2) + (term3 - term4)
  return(likelihood)
}


## Updating slope estimates using Metropolis Hastings MCMC
update_betas<- function(parameters, priors, data){
  
  beta_cand <- rnorm(length(parameters$beta), mean = parameters$beta, sd = priors$beta_prop_sd)
  
  params_top <- list(beta = beta_cand, g = parameters$g, z = parameters$z)
  params_bottom <- list(beta = parameters$beta, g = parameters$g, z = parameters$z)
  
  # Posteriors
  post_top <- loglike(params_top, data) + sum(dnorm(beta_cand, mean = priors$beta_mean, sd = priors$beta_sd, log = TRUE))
  post_bottom <- loglike(params_bottom, data) + sum(dnorm(parameters$beta, mean = priors$beta_mean, sd = priors$beta_sd, log = TRUE))
  
  # Metropolis Hastings Ratio
  log_acceptance_ratio <- post_top - post_bottom
  
  if(log(runif(1)) < log_acceptance_ratio) {
    parameters$beta <- beta_cand
  }
  
  return(parameters)
}

## Updating g(s) using Elliptical Slice Sampling - Calibration Process (used in Source 2 intensity function)
update_g <- function(parameters, priors, data){
  
  # Choosing ellipse (nu) from prior (g)
  nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$g)), Sigma = parameters$tau_2 * exp(-data$dists/4)))
  
  # Log likelihood threshold (finding log(y))
  
  u <- runif(1, min = 0, max = 1)
  
  log_y <- loglike(parameters, data) + log(u)
  
  # Draw an initial proposal for theta
  
  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  
  repeat {
    # Calculate g'
    g_prime <- as.vector((parameters$g - parameters$alpha) *cos(theta) + nu*sin(theta))
    
    params_prime <- parameters
    params_prime$g <- g_prime + parameters$alpha
    
    # Shrinking bracket
    if(loglike(params_prime, data) > log_y){
      parameters$g <- g_prime + parameters$alpha
      return(parameters)
    } else {
      if(theta < 0) {
        theta_min <- theta
      } else {
        theta_max <- theta
      }
      theta <- runif(1, theta_min, theta_max)
    }
  }
}

## Updating z(s) using Elliptical Slice Sampling - Latent Gaussian Random Field
update_z <- function(parameters, priors, data){
  
  # Choosing ellipse (nu) from prior (g)
  nu <- as.vector(MASS::mvrnorm(n = 1, mu = rep(0, length(parameters$z)), Sigma = parameters$sigma_2 * exp(-data$dists/parameters$phi)))
  
  # Log likelihood threshold (finding log(y))
  
  u <- runif(1, min = 0, max = 1)
  
  log_y <- loglike(parameters, data) + log(u)
  
  # Draw an initial proposal for theta
  
  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta - 2*pi
  theta_max <- theta
  
  repeat {
    # Calculate z'
    z_prime <- as.vector(parameters$z*cos(theta) + nu*sin(theta))
    
    params_prime <- parameters
    params_prime$z <- z_prime
    
    # Shrinking bracket
    if(loglike(params_prime, data) > log_y){
      parameters$z <- z_prime
      return(parameters)
    } else {
      if(theta < 0) {
        theta_min <- theta
      } else {
        theta_max <- theta
      }
      theta <- runif(1, theta_min, theta_max)
    }
  }
}

## Updating sigma_2 using Gibbs Sampling - Variance parameter for latent Gaussian random field z(s)
update_sigma_2 <- function(parameters, priors, data){
  n <- length(parameters$z)
  
  alpha_post <- priors$a_0_sigma + n/2
  beta_post  <- priors$b_0_sigma  + 0.5 * t(parameters$z) %*% parameters$Sigma.Inv %*% parameters$z
  
  # Draw samples from Inverse-Gamma
  parameters$sigma_2 <- 1 / rgamma(1, shape = alpha_post, rate = beta_post)
  
  return(parameters)
  
}

## Updating tau_2 using Gibbs Sampling - Variance parameter for calibration process g(s)
update_tau_2 <- function(parameters, priors, data){
  n <- length(parameters$g)
  
  alpha_post <- priors$a_0_tau + n/2
  beta_post  <- priors$b_0_tau  + 0.5 * t(parameters$g - parameters$alpha) %*% parameters$Sigma.Inv %*% (parameters$g - parameters$alpha)
  
  # Draw samples from Inverse-Gamma
  parameters$tau_2 <- 1 / rgamma(1, shape = alpha_post, rate = beta_post)
  
  return(parameters)
  
}

## Updating alpha using Gibbs Sampling - Mean parameter for calibration process g(s)
update_alpha <- function(parameters, priors, data){
  n <- length(parameters$g)
  
  denom <- (sum(parameters$Sigma.Inv) / parameters$tau_2) + (1 / 100)
  mu_post <- (( sum(parameters$Sigma.Inv %*% parameters$g)/ parameters$tau_2 )) / denom
  sigma_post <- sqrt(1 / denom)
  
  parameters$alpha <- rnorm(1, mean = mu_post, sd = sigma_post)
  
  return(parameters)
}

## Wrapper function
driver <- function(parameters, priors, data, iters){
  out=list()
  out$params=parameters
  out$priors=priors
  out$iters=iters
  
  #Posterior containers
  out$beta=matrix(NA,nrow = length(parameters$beta),ncol = iters)
  out$g=matrix(NA, nrow = length(parameters$g), ncol = iters)
  out$z=matrix(NA, nrow = length(parameters$z), ncol = iters)
  out$sigma_2=matrix(NA, nrow = 1, ncol = iters)
  out$tau_2=matrix(NA, nrow = 1, ncol = iters)
  out$alpha=matrix(NA, nrow = 1, ncol = iters)
  
  #Calling each update function
  for(k in 1:iters){
    parameters <- update_betas(parameters, priors, data)
    out$beta[,k]=parameters$beta
    
    parameters <- update_g(parameters, priors, data)
    out$g[,k] <- parameters$g
    
    parameters <- update_z(parameters, priors, data)
    out$z[,k] <- parameters$z
    
    parameters <- update_sigma_2(parameters, priors, data)
    out$sigma_2[,k] <- parameters$sigma_2
    
    parameters <- update_alpha(parameters, priors, data)
    out$alpha[,k] <- parameters$alpha
    
    parameters <- update_tau_2(parameters, priors, data)
    out$tau_2[,k] <- parameters$tau_2
  }
  return(out)
}

# Running Simulation -----------------------------------------------------------------

data <- list(X_grid = X_grid,
             X_1 = X_1,
             X_2 = X_2,
             grid_res = 10,
             cell_area  = (diff(win$xrange) / grid_res) * (diff(win$yrange) / grid_res),
             nn_index_1 = nn_X_1$which,
             nn_index_2 = nn_X_2$which,
             win = win,
             cell_size = cell_size,
             x_seq = x_seq,
             y_seq = y_seq,
             coords = as.matrix(expand.grid(y_seq, x_seq)),
             dists = as.matrix(dist(coords)))

parameters <- list(beta = c(0,0),
                   g = g,
                   z = z,
                   sigma_2 = sigma_2,
                   alpha = alpha,
                   tau_2 = tau_2, 
                   phi = 4)

parameters$Sigma <- exp(-data$dists/ parameters$phi)
parameters$Sigma.Inv <- solve(parameters$Sigma)

priors <- list(beta_mean = c(0,0),
               beta_sd = c(10,10),
               beta_prop_sd = c(0.075, 0.075),
               z_mean = 0,
               a_0_sigma = 3,
               b_0_sigma = 1,
               a_0_tau = 2,
               b_0_tau = 1
)

# save(data, file = "sim_data.RData")
# save(parameters, file = "sim_parameters.RData")
# save(priors, file = "sim_priors.Rdata")

iters <- 10000

burnin <- 3000

sim <- driver(parameters, priors, data, iters) 

# save(sim, file = "sim.RData")

# Trace Plots -------------------------------------------------------------------

plot(sim$beta[1,], type = "l", main = "Beta 1 Trace Plot")
plot(sim$beta[2,], type = "l", main = "Beta 2 Trace Plot")
plot(sim$z[1,], type = "l", main = "z trace plot")
plot(sim$g[1,], type = "l", main = "g trace plot")
plot(sim$sigma_2[1,], type = "l", main = "sigma_2 trace plot")
plot(sim$tau_2[1,], type = "l", main = "tau_2 trace plot")
plot(sim$alpha[1,], type = "l", main = "alpha trace plot")

# Running multiple simulations to calculate empirical coverage ---------------------

n_sims <- 100
n_sims_df <- data.frame()

S <- 0.2 * exp(-dists/4)

mu <- rep(0, n)

covariate <- as.vector(rcpp_rmvnorm(1,S,mu))

for (i in 1:n_sims){
  cat("Running simulation", i, "of", n_sims, "\n")
  
  win <- owin(xrange = c(0, 10), yrange = c(0, 10))
  
  grid_res <- 10
  
  cell_size <- diff(win$xrange) / grid_res
  
  x_seq <- seq(win$xrange[1] + cell_size/2,
               win$xrange[2] - cell_size/2,
               by = cell_size)
  y_seq <- seq(win$yrange[1] + cell_size/2,
               win$yrange[2] - cell_size/2,
               by = cell_size)
  
  coords <- as.matrix(expand.grid(y_seq, x_seq))
  grid_coords <- expand.grid(x = y_seq, y = x_seq)
  
  dists <- as.matrix(dist(coords))
  n <- nrow(coords)
  
  cov_field <- matrix(covariate - mean(covariate),
                      nrow = grid_res,
                      ncol = grid_res,
                      byrow = TRUE)
  
  cov_field_ppp <- ppp(x = grid_coords$x,
                       y = grid_coords$y,
                       window = win,
                       marks = covariate)
  
  #Simulate Gaussian random field
  
  sigma_2 <- 0.10
  
  S_z <-  sigma_2 * exp(-dists/4)
  
  z <- as.vector(rcpp_rmvnorm(1,S_z,mu))
  
  z_mat <- matrix(z, nrow = length(x_seq), ncol = length(y_seq))
  
  z_ppp <- ppp(x = grid_coords$x,
               y = grid_coords$y,
               window = win,
               marks = z)
  
  # Simulate LGCP
  b_0 <- 1.75
  b_1 <- 3
  
  lambda <- exp(b_0 + b_1*(cov_field) + z)
  lambda_im <- im(lambda, xcol = x_seq, yrow = y_seq)
  
  lgcp_sim <- rpoispp(lambda_im)
  lgcp_sim
  
  plot(lgcp_sim, main = paste("lgcp_sim", i))
  
  # Discretize using spatstat
  
  lgcp_discretize <- pixellate(lgcp_sim, eps = 1)
  
  # Source 1  -----------
  nrow <- length(lgcp_discretize$yrow)
  ncol <- length(lgcp_discretize$xcol)
  
  x_min_subwindow1 <- 0
  x_max_subwindow1 <- 5
  y_min_subwindow1 <- 5
  y_max_subwindow1 <- 10
  
  x_min_subwindow2 <- 5
  x_max_subwindow2 <- 10
  y_min_subwindow2 <- 0
  y_max_subwindow2 <- 5
  
  sub_window1 <- owin(xrange = c(x_min_subwindow1, x_max_subwindow1), yrange = c(y_min_subwindow1, y_max_subwindow1))
  
  sub_window2 <- owin(xrange = c(x_min_subwindow2, x_max_subwindow2), yrange = c(y_min_subwindow2, y_max_subwindow2))
  
  lambda_1 <- lgcp_discretize
  
  lgcp_1_sub1 <- rpoispp(lambda_1[sub_window1])
  lgcp_1_sub2 <- rpoispp(lambda_1[sub_window2])
  
  lgcp_1 <- superimpose(lgcp_1_sub1, lgcp_1_sub2, W = owin(c(0,10), c(0,10)))
  
  # Source 2 --------
  
  tau_2 <- 0.4
  S_g <- tau_2 * exp(-dists/4)
  alpha <- -0.2
  
  #g <- rnorm(nrow * ncol, alpha, tau_2)
  g <- as.vector(rcpp_rmvnorm(1,S_g,alpha))
  
  exp_g <- exp(g)
  
  lambda_2 <- lgcp_discretize * exp_g
  
  lgcp_2 <- rpoispp(lambda_2)
  
  # Data frame creation -----
  
  ## X grid
  X_grid <- as.data.frame(coords)
  colnames(X_grid) <- c("x", "y")
  X_grid$covariate <- covariate
  
  ## Source 1 (small var)
  X_1 <- as.data.frame(lgcp_1)
  nn_X_1 <- nncross(lgcp_1, cov_field_ppp)
  X_1$covariate <- cov_field_ppp$marks[nn_X_1$which]
  
  # making a mask for X-grid to be used with source 1
  inside_sub1 <- with(X_grid,
                      x >= x_min_subwindow1 & x <= x_max_subwindow1 &
                        y >= y_min_subwindow1 & y <= y_max_subwindow1)
  
  inside_sub2 <- with(X_grid,
                      x >= x_min_subwindow2 & x <= x_max_subwindow2 &
                        y >= y_min_subwindow2 & y <= y_max_subwindow2)
  
  X_grid$mask_source1 <- inside_sub1 | inside_sub2
  
  ## Source 2 (large var)
  X_2 <- as.data.frame(lgcp_2)
  nn_X_2 <- nncross(lgcp_2, cov_field_ppp)
  X_2$covariate <- cov_field_ppp$marks[nn_X_2$which]
  
  data <- list(X_grid = X_grid,
               X_1 = X_1,
               X_2 = X_2,
               grid_res = 10,
               cell_area  = (diff(win$xrange) / grid_res) * (diff(win$yrange) / grid_res),
               nn_index_1 = nn_X_1$which,
               nn_index_2 = nn_X_2$which,
               win = win,
               cell_size = cell_size,
               x_seq = x_seq,
               y_seq = y_seq,
               coords = as.matrix(expand.grid(y_seq, x_seq)),
               dists = as.matrix(dist(coords)))
  
  parameters <- list(beta = c(0,0),
                     g = g,
                     z = z,
                     sigma_2 = sigma_2,
                     alpha = alpha,
                     tau_2 = tau_2, 
                     phi = 4)
  
  parameters$Sigma <- exp(-data$dists/ parameters$phi)
  parameters$Sigma.Inv <- solve(parameters$Sigma)
  
  priors <- list(beta_mean = c(0,0),
                 beta_sd = c(10,10),
                 beta_prop_sd = c(0.075, 0.075),
                 z_mean = 0,
                 a_0_sigma = 3,
                 b_0_sigma = 1,
                 a_0_tau = 2,
                 b_0_tau = 1
  )
  
  iters <- 10000
  
  burnin <- 3000
  
  # Run simulation ----------
  
  sim <- driver(parameters, priors, data, iters)
  
  beta_post <- sim$beta[, (burnin+1):iters]
  alpha_post <- sim$alpha[, (burnin+1):iters]
  sigma_2_post <- sim$sigma_2[, (burnin+1):iters]
  tau_2_post <- sim$tau_2[, (burnin+1):iters]
  g_post <- sim$g[,(burnin+1):iters]
  z_post <- sim$z[,(burnin+1):iters]
  
  # Calculating posterior
  posterior_lambda <- matrix(NA, nrow = nrow(X_grid), ncol = (iters-burnin))
  
  for(m in 1:(iters-burnin)){
    beta_m <- beta_post[,m]
    
    log_lambda_m <- beta_m[1] + beta_m[2]*covariate + z_post[,m]
    
    posterior_lambda[, m] <- exp(log_lambda_m)
  }
  
  # Posterior mean intensity
  lambda_mean <- rowMeans(posterior_lambda, na.rm = TRUE)
  
  actual_count <- sum(lambda * data$cell_area, na.rm = TRUE)
  
  n_draws <- ncol(posterior_lambda)
  expected_total_per_draw_fused <- numeric(n_draws)
  
  for (m in 1:n_draws) {
    expected_total_per_draw_fused[m] <- sum(posterior_lambda[, m] * data$cell_area, na.rm = TRUE)
  }
  
  expected_mean_fused <- mean(expected_total_per_draw_fused)
  expected_sd_fused   <- sd(expected_total_per_draw_fused)
  
  # Store just summary stats
  n_sims_df <- rbind(n_sims_df, data.frame(
    sim = i,
    sim_n = lgcp_sim$n,
    
    # beta parameters
    beta0_mean = mean(beta_post[1,]), 
    beta0_sd = sd(beta_post[1,]),
    beta0_lower = quantile(beta_post[1,], 0.025),
    beta0_upper = quantile(beta_post[1,], 0.975),
    
    beta1_mean = mean(beta_post[2,]), 
    beta1_sd = sd(beta_post[2,]),
    beta1_lower = quantile(beta_post[2,], 0.025),
    beta1_upper = quantile(beta_post[2,], 0.975),
    
    # hyperparameters
    sigma2_mean = mean(sigma_2_post), 
    sigma2_sd = sd(sigma_2_post),
    sigma2_lower = quantile(sigma_2_post, 0.025),
    sigma2_upper = quantile(sigma_2_post, 0.975),
    
    tau2_mean = mean(tau_2_post), 
    tau2_sd = sd(tau_2_post),
    tau2_lower = quantile(tau_2_post, 0.025),
    tau2_upper = quantile(tau_2_post, 0.975),
    
    alpha_mean = mean(alpha_post), 
    alpha_sd = sd(alpha_post),
    alpha_lower = quantile(alpha_post, 0.025),
    alpha_upper = quantile(alpha_post, 0.975),
    
    g_mean = mean(g_post),
    g_sd = sd(g_post),
    g_lower = quantile(g_post, 0.025),
    g_upper = quantile(g_post, 0.975),
    
    z_mean = mean(z_post),
    z_sd = sd(z_post),
    z_lower = quantile(z_post, 0.025),
    z_upper = quantile(z_post, 0.975),
    
    #expected counts
    count_mean = expected_mean_fused,
    count_lower = quantile(expected_total_per_draw_fused, 0.025),
    count_upper = quantile(expected_total_per_draw_fused, 0.975),
    actual_count = actual_count
    
  ))
}

# save(n_sims_df, file = "n_sims_df.Rdata")

