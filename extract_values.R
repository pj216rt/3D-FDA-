# source("data_gen.R")
source("STAN_control.R")

#want a function to extract the values of various values
fixed.effects.analysis <- function(stan_object, ci.level=0.95){
  #specifcying the alpha for the given CI level
  alpha <- (1 - ci.level) / 2
  BETA.results <- stan_object$summary(
    variables = "BETA",
    mean,
    ~quantile(.x, probs = c(alpha, 1 - alpha))
  )
  
  #rename columns to clarify
  colnames(BETA.results) <- c("variable", "mean", "CI_lower", "CI_upper")
  BETA.results <- BETA.results %>%
    mutate(
      parameter = as.integer(gsub("BETA\\[([0-9]+),([0-9]+)\\]", "\\1", variable)),  
      grid_point = as.integer(gsub("BETA\\[([0-9]+),([0-9]+)\\]", "\\2", variable)) 
    )
  
  n_predictors <- n_distinct(BETA.results$parameter)
  n_dims <- 3 #three dimensional
  n_per_dim <- n_distinct(BETA.results$grid_point)/n_dims
  print(n_per_dim)
  
  beta_summary <- BETA.results %>%
    mutate(
      predictor = parameter,
      dim = (grid_point - 1) %/% n_per_dim + 1,
      position = (grid_point - 1) %% n_per_dim + 1) %>%
    dplyr::select(predictor, dim, position, mean, CI_lower, CI_upper)
  
  return(beta_summary)
}

summarized_beta_values <- fixed.effects.analysis(stan_object = results)

#function to plot this data
plot.summary <- function(beta.summarized){
  stan.plot <- ggplot(data = beta.summarized, 
                      aes(x = position, y = mean, color = as.factor(predictor), group = predictor)) + 
    geom_point(size = 0.3) + 
    geom_line() + 
    facet_wrap(~ dim, scales = "free", strip.position = "top") + 
    guides(col = guide_legend(title = "Predictors")) + 
    theme_classic() + 
    labs(title = "Posterior Means STAN", 
         x = "Time Grid Point", 
         y = "Value") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2)  # Ensure correct CI columns
  
  #should save this
  file_name = "recovered_means_plots.jpg"
  ggsave(stan.plot, file=file_name, path = "plots/post_plots/")
  
  return(stan.plot)
}

summarized.plots <- plot.summary(beta.summarized = summarized_beta_values)

#function to compare the true beta curves with the recovered beta curves
# plot.compare <- function(path1, path2){
#   #load plots
#   img1 <- readJPEG(path1)
#   img2 <- readJPEG(path2)
#   
#   combined <- grid.arrange(rasterGrob(img1), rasterGrob(img2), ncol = 2)
#   print(combined)
# }

# comparison.plots <- plot.compare(path1="plots/beta_plots/beta_plot.jpg",
#                      path2="plots/post_plots/recovered_means_plots.jpg")

#function compute coverage of coefficients
#can use the summarized_beta_values and see how many times for each
#If the true value of the curve at a given grid point is true_value, 
#and the CI is [CI_lower, CI_upper], then:
#Coverage = 1 if true_value is within [CI_lower, CI_upper]
#Coverage = 0 if true_value is outside this range.
#Then you can compute the proportion of grid points where the true value lies within the CI

#coverage calculation:
compute.coverage <- function(true_values, summary.dat, grid.size){
  #pivot true values
  true_values_long <- true_values %>%
    as.data.frame() %>%
    mutate(predictor = 1:nrow(true_values)) %>% 
    pivot_longer(cols = starts_with("V"), names_to = "position", 
                 values_to = "true_value") %>%
    mutate(position = as.integer(gsub("V", "", position)),  
           dim = (position - 1) %/% grid.size + 1) %>%
    group_by(dim, predictor, position)
  
  coverage_data <- summary.dat %>%
    inner_join(true_values_long, by = c("predictor", "position", "dim")) %>%
    mutate(
      #Check if the true value lies within the CI bounds
      coverage = ifelse(true_value >= CI_lower & true_value <= CI_upper, 1, 0)
    )
  
  #Compute the coverage rate for each predictor
  coverage_summary <- coverage_data %>%
    group_by(predictor) %>%
    summarise(coverage_rate = mean(coverage))
  
  return(coverage_summary)
}

test <- compute.coverage(true_values = as.data.frame(beta.coefficients), 
                         summary.dat = summarized_beta_values, grid.size = 25)

#need function to extract MSE of Y
compute.mse.y <- function(results){
  #get MSE_Y samples
  MSE_Y_samples <- results$draws(variables = "MSE_Y")
  
  #compute mean and 95% CI
  MSE_Y_summary <- summarise_draws(MSE_Y_samples, 
                                   mean,
                                   ~quantile(.x, probs = c(0.025, 0.975)))
  colnames(MSE_Y_summary) <- c("variable", "mean", "CI_lower", "CI_upper")
  
  return(MSE_Y_summary)
}

#function to compute integrated mean squared error
compute.IMSE <- function(summary.dat, true_values, grid.size) {
  #pivot true values
  true_values_long <- true_values %>%
    as.data.frame() %>%
    mutate(predictor = 1:nrow(true_values)) %>% 
    pivot_longer(cols = starts_with("V"), names_to = "position", 
                 values_to = "true_value") %>%
    mutate(position = as.integer(gsub("V", "", position)),  
           dim = (position - 1) %/% grid.size + 1) %>%
    group_by(dim, predictor, position) %>% group_by(dim, position, predictor)

  # #join the data to compute squared difference, then pivot wider so that we
  # #have a p*3D matrix of the squared differences
  MSE.dat <- summary.dat %>%
    inner_join(true_values_long, by = c("predictor", "position", "dim")) %>% 
    mutate(squared_diff = (mean - true_value)^2)

  #Calculate the IMSe for each predictor
  IMSE_vals <- MSE.dat %>% group_by(predictor) %>%
    summarise(IMSE = mean(squared_diff))


  return(IMSE_vals)
}

test1 <- compute.IMSE(summary.dat = summarized_beta_values, true_values = beta.coefficients,
                      grid.size = ncol(beta.coefficients)/3)
head(test1)
