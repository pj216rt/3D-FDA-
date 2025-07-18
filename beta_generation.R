library(tidyverse)
library(mvtnorm)
library(spam)  #not entirely sure what the spam call in the sampler does
library(splines)
library(MASS)
library(MCMCpack) #for the inverse wishart distribution
require(gridExtra)
library(bayesplot)
library(grid)

#function to generate smooth functional coefficients using B-splines
beta.gen.fun.tri <- function(grid.size, knots, parameters, deg){
  x <- seq(0, grid.size, length.out = grid.size)
  
  #b spline basis matrix
  bs_basis <- bs(x, degree = deg)  
  
  beta_list <- vector("list", 3)
  
  #going to loop over each dimension
  for (dim in 1:3) {
    random_coeffs_matrix <- matrix(rnorm(ncol(bs_basis)*parameters), 
                                   nrow = ncol(bs_basis), ncol = parameters)
    
    #compute the beta matrix for dimension dim
    beta_matrix <- bs_basis %*% random_coeffs_matrix
    
    #store this in the list
    beta_list[[dim]] <- t(beta_matrix)
  }
  
  #combine all three matrices column wise
  coefs <- do.call(cbind, beta_list)
  
  return(coefs)
}

beta.coefficients <- beta.gen.fun.tri(grid.size = 25, knots = 8, parameters = 3, deg=3)

#need a function to plot these generated coefficients
plot.beta.coef <- function(beta.matrix){
  #compute grid size in 1 dimension
  grid.size <- ncol(beta.matrix)/3
  
  #parameters
  params <- nrow(beta.matrix)
  
  df_list <- list()
  for (i in 1:3){
    #need to extract "blocks" of columns from the beta coefficient matrix
    #selecting the starting and stopping columns to extract for each dim
    start.col <- (i - 1)*grid.size + 1
    end.col <- i*grid.size
    
    #extracting elements of the beta coefficients and naming them
    #adding an identifier for parameter and dimension
    df_temp <- as.data.frame(beta.matrix[, start.col:end.col])
    colnames(df_temp) <- paste0("Grid_", 1:grid.size)
    df_temp$Parameter <- factor(1:params)
    df_temp$Dimension <- factor(i)

    df_list[[i]] <- df_temp
  }
  
  #combine the elements of this list
  df_long <- do.call(rbind, df_list)
  
  #need to convert to long format for ggplot
  #need to convert the grid to numeric
  df_long <- df_long %>%
    pivot_longer(cols = starts_with("Grid_"), 
                 names_to = "Grid", 
                 values_to = "Beta_Value") %>%
    mutate(Grid = as.numeric(sub("Grid_", "", Grid)))
  
  #now use ggplot to display
  #labels for the facet labeller
  df_long$Dimension <- factor(df_long$Dimension, levels = c(1, 2, 3), 
                              labels = c("Dimension X", "Dimension Y", "Dimension Z"))
  
  beta.plot <- ggplot(df_long, aes(x = Grid, y = Beta_Value, color = Parameter, group = Parameter)) +
    geom_line() +
    facet_wrap(~ Dimension, ncol = 1) + 
    labs(title = "Functional Coefficients Across Dimensions",
         x = "Grid Points",
         y = "Coefficient Value") + 
    theme_bw()
  
  print(beta.plot)
  
  #need to save this plot
  file_name = "beta_plot.jpg"
  ggsave(beta.plot, file=file_name, path = "plots/beta_plots/")
  
}

beta.coef.plots <- plot.beta.coef(beta.matrix = beta.coefficients)
