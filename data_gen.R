#load in the beta generation stuff
source("beta_generation.R")

#generating 3 dimenaional data, assuming balanced design
gen.3d.data <- function(num.subj = 50, num.visits = 5,
                        beta.vals, var.x = 1, var.z=1){
  #make sure you supply beta vals
  if (missing(beta.vals)) {
    message("Beta values not supplied.")
  }
  
  grid.size.3d <- ncol(beta.vals)
  n.params <-nrow(beta.vals)
  
  #total observations
  tot.obs <- num.subj*num.visits
  
  #generate the X design matrix
  X <- matrix(rnorm(n=(num.subj*n.params), sd = sqrt(var.x)), 
              nrow = num.subj, ncol = n.params) %>%
    as.data.frame() %>%
    mutate(subj = row_number())
  X$subj = as.factor(X$subj)
  X1 <- X[rep(seq_len(nrow(X)), each = num.visits), ]
  X.des <- model.matrix(~. -1 -subj, data = X1)
  
  #generate fixed effects
  fixef = as.matrix(X.des) %*% as.matrix(beta.vals)
  
  #random effects
  Z.des = model.matrix( ~ 0 + subj + (-1):subj, data = X1)
  subj.ranef <- matrix(rnorm(n=(num.subj*grid.size.3d), sd=sqrt(var.z)), 
                       nrow = num.subj, ncol = grid.size.3d)
  ranef <- Z.des %*% subj.ranef
  
  #level 1 residuals
  eps <- rnorm(n = tot.obs)
  
  Yij.true <- fixef + ranef
  Yij.obs <- fixef + ranef + eps
  
  #return "observed" data and other stuff in a list
  
  output.list <- list()
  output.list$obs <- Yij.obs
  output.list$true <- Yij.true
  output.list$x_design <- X.des
  output.list$z_design <- Z.des
  output.list$raw.data <- X1
  return(output.list)
  
}

generated.data <- gen.3d.data(beta.vals=beta.coefficients)

#function to plot generated data
plot.3d.generated.dat <- function(num.subj=50, num.visits=5, observed.3d.dat){
  #get the grid size of one dimension
  grid.size <- ncol(observed.3d.dat$obs) / 3
  
  #get dataframe of observed data
  df_wide <- as.data.frame(observed.3d.dat$obs)
  df_wide <- df_wide %>%
    mutate(Subject = rep(1:num.subj, each = num.visits),
           Visit = rep(1:num.visits, num.subj))
  
  #ppivot this
  df_long <- df_wide %>%
    pivot_longer(cols = matches("^V[0-9]+$"),  
                 names_to = "GridIndex", 
                 values_to = "FunctionValue") %>%
    mutate(GridIndex = as.numeric(str_remove(GridIndex, "V")),
           Dimension = case_when(
             GridIndex <= grid.size ~ "Dimension X",
             GridIndex <= 2 * grid.size ~ "Dimension Y",
             GridIndex <= 3 * grid.size ~ "Dimension Z",
             TRUE ~ "Unknown"),
           GridPoint = (GridIndex - 1) %% grid.size + 1)
  
  #plot this data
  obs.dat.plots <- ggplot(df_long, aes(x = GridPoint, y = FunctionValue, 
                                       group = interaction(Subject, Visit), 
                                       color = as.factor(Subject))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~Dimension, scales = "free_y") +  # Separate plots for each dimension
    labs(title = "Simulated 3D Functional Data", 
         x = "Grid Point", 
         y = "Function Value", 
         color = "Subject") +
    theme_minimal() +
    theme(legend.position = "none")
  print(obs.dat.plots)
  
  #need to save this plot
  file_name = "obs_data_plot.jpg"
  ggsave(obs.dat.plots, file=file_name, path = "plots/data_plots/")
}

observed.data.plots <- plot.3d.generated.dat(observed.3d.dat = generated.data)

