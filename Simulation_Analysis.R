# Required Packages -------------------------------------------------------------
library(spatstat)
library(fields)
library(FastGP)
library(MASS)
library(ggplot2)
library(coda)

# Simulation data from Data_Fusion_Simulation_Study.R ---------------------------
load("sim_full.RData") # Simulation using model fit to both Source 1 and Source 2 data
load("sim_source1.RData") # Simulation using model fit to only Source 1 data
load("sim_source2.RData") # Simulation using model fit to only Source 2 data
load("sim_data.RData") # List of data used in simulation
load("sim_parameters.RData") # List of parameters used in simulation
load("sim_priors.Rdata") # List of priors used in simulation
load("sim_covariate.Rdata") # vector of simulated covariate values
load("sim_lambda.RData") # Data frame of intensity values for each pixel within domain
load("n_sims_df.RData") # 100 iterations of simulation study with 95% credible intervals for each parameter

# Domain of simulation study ----------------------------------------------------
grid_res <- 10

cell_size <- diff(win$xrange) / grid_res

x_seq <- seq(win$xrange[1] + cell_size/2,
             win$xrange[2] - cell_size/2,
             by = cell_size)
y_seq <- seq(win$yrange[1] + cell_size/2,
             win$yrange[2] - cell_size/2,
             by = cell_size)
grid_coords <- expand.grid(x = y_seq, y = x_seq)

burnin <- 3000
iters <- 10000

# Mean, standard deviation, and 95% credible intervals for each parameter --------

## Simulation using model fit to both Source 1 and Source 2 
beta_post <- sim$beta[, (burnin+1):iters]
alpha_post <- sim$alpha[,(burnin+1):iters]
sigma_2_post <- sim$sigma_2[,(burnin+1):iters]
tau_2_post <- sim$tau_2[,(burnin+1):iters]
g_post <- sim$g[,(burnin+1):iters]
z_post <- sim$z[,(burnin+1):iters]

apply(beta_post, 1, mean)
apply(beta_post, 1, sd)
quantile(beta_post[1,], c(0.025, 0.975))
quantile(beta_post[2,], c(0.025, 0.975))

mean(sigma_2_post)
sd(sigma_2_post)
quantile(sigma_2_post, c(.025,0.975))

mean(tau_2_post)
sd(tau_2_post)
quantile(tau_2_post, c(.025,0.975))

mean(alpha_post)
sd(alpha_post)
quantile(alpha_post, c(.025,0.975))

mean(g_post)
sd(g_post)
quantile(g_post, c(.025,0.975))

mean(z_post)
sd(z_post)
quantile(z_post, c(.025,0.975))

## Simulation using model fit to only Source 1 data
beta_post1 <- sim_source1$beta[, (burnin+1):iters]
alpha_post1 <- sim_source1$alpha[,(burnin+1):iters]
sigma_2_post1 <- sim_source1$sigma_2[,(burnin+1):iters]
tau_2_post1 <- sim_source1$tau_2[,(burnin+1):iters]
g_post1 <- sim_source1$g[,(burnin+1):iters]
z_post1 <- sim_source1$z[,(burnin+1):iters]

apply(beta_post1, 1, mean)
apply(beta_post1, 1, sd)
quantile(beta_post1[1,], c(0.025, 0.975))
quantile(beta_post1[2,], c(0.025, 0.975))

mean(sigma_2_post1)
sd(sigma_2_post1)
quantile(sigma_2_post1, c(.025,0.975))

mean(z_post1)
sd(z_post1)
quantile(z_post1, c(.025,0.975))


## Simulation using model fit to only Source 2 data
beta_post2 <- sim_source2$beta[, (burnin+1):iters]
alpha_post2 <- sim_source2$alpha[,(burnin+1):iters]
sigma_2_post2 <- sim_source2$sigma_2[,(burnin+1):iters]
tau_2_post2 <- sim_source2$tau_2[,(burnin+1):iters]
g_post2 <- sim_source2$g[,(burnin+1):iters]
z_post2 <- sim_source2$z[,(burnin+1):iters]

apply(beta_post2, 1, mean)
apply(beta_post2, 1, sd)
quantile(beta_post2[1,], c(0.025, 0.975))
quantile(beta_post2[2,], c(0.025, 0.975))

mean(sigma_2_post2)
sd(sigma_2_post2)
quantile(sigma_2_post2, c(.025,0.975))

mean(z_post2)
sd(z_post2)
quantile(z_post2, c(.025,0.975))

# Posterior Intensity Surfaces --------------------------------------------------

## Posterior intensity surface for fusion model
posterior_lambda <- matrix(NA, nrow = nrow(data$X_grid), ncol = (iters-burnin))

for(m in 1:(iters-burnin)){
  beta_m <- beta_post[,m]
  
  log_lambda_m <- beta_m[1] + beta_m[2]*covariate + z_post[,m]
  
  posterior_lambda[, m] <- exp(log_lambda_m)
}

lambda_mean <- rowMeans(posterior_lambda, na.rm = TRUE)

lambda_mean_mat <- matrix(lambda_mean, 
                          nrow = data$grid_res, 
                          ncol = data$grid_res, 
                          byrow = TRUE)

## Posterior intensity surface for model fit to only source 1 data
posterior_lambda_source1 <- matrix(NA, nrow = nrow(data$X_grid), ncol = (iters-burnin))

for(m in 1:(iters-burnin)){
  beta_m1 <- beta_post1[,m]
  
  log_lambda_m1 <- beta_m1[1] + beta_m1[2]*covariate + z_post1[,m]
  
  posterior_lambda_source1[, m] <- exp(log_lambda_m1)
}

lambda_mean1 <- rowMeans(posterior_lambda_source1, na.rm = TRUE)

lambda_mean_mat1 <- matrix(lambda_mean1, 
                           nrow = data$grid_res, 
                           ncol = data$grid_res, 
                           byrow = TRUE)

## Posterior intensity surface for model fit to only source 2 data
posterior_lambda_source2 <- matrix(NA, nrow(data$X_grid), ncol = (iters-burnin))

for(m in 1:(iters-burnin)){
  beta_m2 <- beta_post2[,m]
  
  log_lambda_m2 <- beta_m2[1] + beta_m2[2]*covariate + z_post2[,m]
  
  posterior_lambda_source2[, m] <- exp(log_lambda_m2)
}

lambda_mean2 <- rowMeans(posterior_lambda_source2, na.rm = TRUE)

lambda_mean_mat2 <- matrix(lambda_mean2, 
                           nrow = data$grid_res, 
                           ncol = data$grid_res, 
                           byrow = TRUE)

## Defining common scale for all posterior intensity surfaces
all_log_values <- c(
  log(lambda),
  log(lambda_mean),
  log(lambda_mean1),
  log(lambda_mean2)
)

zlim <- range(all_log_values, na.rm = TRUE)

## Plotting true posterior intensity surfaces

# pdf(file = "true_intensity.pdf", width = 5, height = 5)
quilt.plot(grid_coords$y,
           grid_coords$x,
           log(lambda), 
           nx = 10,
           ny = 10,
           col = viridis(100),
           zlim = zlim,
           main = "True Intensity Surface")
# dev.off()

# pdf(file = "fused_intensity_sim.pdf", width = 5, height = 5)
quilt.plot(grid_coords$x,
           grid_coords$y,
           log(lambda_mean),
           nx = 10, 
           ny = 10, 
           col = viridis(100),
           zlim = zlim,
           main = "Fusion Model Posterior Intensity Surface")
# dev.off()

# pdf(file = "source1_intensity_sim.pdf", width = 5, height = 5)
quilt.plot(grid_coords$x,
           grid_coords$y,
           log(lambda_mean1),
           nx = 10, 
           ny = 10, 
           col = viridis(100),
           zlim = zlim,
           main = "Source 1 Posterior Intensity Surface")
# dev.off()

# pdf(file = "source2_intensity_sim.pdf", width = 5, height = 5)
quilt.plot(grid_coords$x,
           grid_coords$y,
           log(lambda_mean2),
           nx = 10, 
           ny = 10, 
           col = viridis(100),
           zlim = zlim,
           main = "Source 2 Posterior Intensity Surface")
# dev.off()

## Plotting differences in true intensity surface and posterior intensity surfaces

z_lim_diff <- range(c(log(lambda_mean) - t(log(lambda)), 
                      log(lambda_mean1) - t(log(lambda)), 
                      log(lambda_mean2) - t(log(lambda))))

# pdf(file = "diff_fusion.pdf", width = 5, height = 5)
quilt.plot(grid_coords$x,
           grid_coords$y,
           log(lambda_mean) - t(log(lambda)) ,
           nx = 10,
           ny = 10,
           zlim = z_lim_diff,
           col = colorRampPalette(colors = c("red", "white", "blue"))(100),
           main = "Difference in True and \n Fusion Model Log-Intensity Surface")
# dev.off()

# pdf(file = "diff_source1.pdf", width = 5, height = 5)
quilt.plot(grid_coords$x,
           grid_coords$y,
           log(lambda_mean1) - t(log(lambda)),
           nx = 10,
           ny = 10, 
           zlim = z_lim_diff,
           col = colorRampPalette(colors = c("red", "white", "blue"))(100),
           main = "Difference in True and \n Source 1 Model Log-Intensity Surface")
# dev.off()

# pdf(file = "diff_source2.pdf", width = 5, height = 5)
quilt.plot(grid_coords$x,
           grid_coords$y,
           log(lambda_mean2) - t(log(lambda)),
           nx = 10,
           ny = 10, 
           zlim = z_lim_diff, 
           col = colorRampPalette(colors = c("red", "white", "blue"))(100), 
           main = "Difference in True and \n Source 2 Log-Intensity Surface")
# dev.off()

# Posterior calibration process surface -----------------------------------------
posterior_g <- matrix(NA, nrow = nrow(data$X_grid), ncol = (iters - burnin))

for (m in 1:(iters - burnin)) {
  log_g_m <-  g_post[, m]
  posterior_g[, m] <- exp(log_g_m)
}

g_mean <- rowMeans(posterior_g)

## Plotting posterior calibration process surface
# pdf(file = "g_intensity_sim.pdf", width = 5, height = 5)
quilt.plot(grid_coords$y, 
           grid_coords$x, 
           log(g_mean),
           nx = 10,
           ny = 10, 
           main = "Posterior Calibration Parameter Surface",
           col = viridis(100))
# dev.off()

# Analysis of 100 simulations for empirical coverage calculations ----------------
## Beta 
covered_beta0 <- mean(n_sims_df$beta0_lower <= 1.5 &
                        n_sims_df$beta0_upper >= 1.5)

covered_beta1 <- mean(n_sims_df$beta1_lower <= 3 &
                        n_sims_df$beta1_upper >= 3)
covered_beta0
covered_beta1

## sigma_2
covered_sigma2 <- mean(n_sims_df$sigma2_lower <= 0.10 &
                         n_sims_df$sigma2_upper >=0.10)
covered_sigma2

#tau_2
covered_tau2 <- mean(n_sims_df$tau2_lower <= 0.4 &
                       n_sims_df$tau2_upper >=0.4)
covered_tau2

#alpha

covered_alpha <- mean(n_sims_df$alpha_lower <= -0.2 &
                        n_sims_df$alpha_upper >= -0.2)
covered_alpha

#covered count

covered_count <- mean(n_sims_df$count_lower <= n_sims_df$actual_count &
                        n_sims_df$count_upper >= n_sims_df$actual_count)

covered_count

# Expected counts -------------------------------------------------------

## True expected count
actual_count <- sum(lambda * data$cell_area, na.rm = TRUE)
actual_count

n_draws <- ncol(posterior_lambda)

## Expected count from posterior fusion model intensity function
expected_total_per_draw_fused <- numeric(n_draws)

for (m in 1:n_draws) {
  expected_total_per_draw_fused[m] <- sum(posterior_lambda[, m] * data$cell_area, na.rm = TRUE)
}

expected_mean_fused <- mean(expected_total_per_draw_fused)
expected_sd_fused   <- sd(expected_total_per_draw_fused)

count_lower <-  quantile(expected_total_per_draw_fused, 0.025)
count_upper <- quantile(expected_total_per_draw_fused, 0.975)

expected_mean_fused
expected_sd_fused 
count_lower
count_upper

## Expected count from posterior source 1 intensity function
expected_total_per_draw_fused1 <- numeric(n_draws)

for (m in 1:n_draws) {
  expected_total_per_draw_fused1[m] <- sum(posterior_lambda_source1[, m] * data$cell_area, na.rm = TRUE)
}

expected_mean_fused1 <- mean(expected_total_per_draw_fused1)
expected_sd_fused1  <- sd(expected_total_per_draw_fused1)

count_lower1 <-  quantile(expected_total_per_draw_fused1, 0.025)
count_upper1 <- quantile(expected_total_per_draw_fused1, 0.975)

expected_mean_fused1
expected_sd_fused1
count_lower1
count_upper1

## Expected count from posterior source 2 intensity function
expected_total_per_draw_fused2 <- numeric(n_draws)

for (m in 1:n_draws) {
  expected_total_per_draw_fused2[m] <- sum(posterior_lambda_source2[, m] * data$cell_area, na.rm = TRUE)
}

expected_mean_fused2 <- mean(expected_total_per_draw_fused2)
expected_sd_fused2  <- sd(expected_total_per_draw_fused2)

count_lower2 <-  quantile(expected_total_per_draw_fused2, 0.025)
count_upper2 <- quantile(expected_total_per_draw_fused2, 0.975)

expected_mean_fused2
expected_sd_fused2
count_lower2
count_upper2

# RMSE calculations ---------------------------------------------------------------

rmse_fusion <- sqrt(1/100 * (sum(lambda - lambda_mean_mat))^2)

rmse_source_1 <- sqrt(1/100 * (sum(lambda - lambda_mean_mat1))^2)

rmse_source_2 <- sqrt(1/100 * (sum(lambda - lambda_mean_mat2))^2)

print(c(rmse_fusion, rmse_source_1, rmse_source_2))

