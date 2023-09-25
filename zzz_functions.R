# Function to plot cleveland dotplots for all continuous variables in a data.frame
outlier_check <- function(x, variables = NULL, dot_color = "steelblue", dot_size = 3) {
  if (is.null(variables)) {
    x_numeric <-
      x %>% 
      select_if(is.numeric) 
    variables <- colnames(x_numeric)
  }
  
  x_numeric <-
    x %>% 
    select_if(is.numeric) %>% 
    select(variables) 
  
  x_numeric$order <- seq_len(nrow(x_numeric))
  
  x_names  <- colnames(x_numeric)
  x_labels <- labels
  x_longer <- pivot_longer(x_numeric, -order)
  
  # Plot facet
  ggplot(x_longer, aes(y = value, x = order)) +
    geom_point(size = dot_size, alpha = .5, fill = dot_color, shape = 21, color = "black") +
    facet_wrap(~ name, scales = "free") +
    coord_flip() +
    labs(y = "Value of the variable",
         x = "Order of the data") +
    theme_classic(base_size = 18) +
    ggtitle("Outliers check with Cleveland dotplot")
} 

# Function to plot histograms fro all continous variables in a data.frame
normality_check <- function(x, variables = NULL, dot_color = "steelblue", title = "Normality check with histograms") {
  if (is.null(variables)) {
    x_numeric <-
      x %>% 
      select_if(is.numeric) 
    variables <- colnames(x_numeric)
  }
  
  x_numeric <-
    x %>% 
    select_if(is.numeric) %>% 
    select(variables) 
  
  x_numeric$order <- seq_len(nrow(x_numeric))
  
  x_names  <- colnames(x_numeric)
  x_labels <- labels
  x_longer <- pivot_longer(x_numeric, -order)
  
  # Plot facet
  ggplot(x_longer, aes(x = value)) +
    geom_histogram(color = "black", fill = dot_color) +
    facet_wrap(~ name, scales = "free") +
    labs(y = "Frequency",
         x = "Value of the variable") +
    theme_classic(base_size = 18) +
    ggtitle(title)
} 

# Function to check the excess of zeroes
zero_check <- function(x, variable, n_simu = 1000) {

  
  y <- pull(x, variable)
  n <- length(y)
  sim_zip <- numeric()
  # observed zero-inflation idex
  p0_obs  <- sum(y == 0)/n
  zip_obs <- 1 + log(p0_obs) / mean(y)
  
  for (i in 1:n_simu){
    simdata <- rpois(lambda =  mean(y), n = nrow(x))
    sim_p0  <- sum(simdata == 0) / length(simdata)
    if(sim_p0 == 0) {
      sim_p0 <- 1/n
    }
    sim_zip[i] <- 1 + log(sim_p0) / mean(simdata)
  }
  
  res <- data.frame(zip = sim_zip, N = seq_len(length(sim_zip)))
  p_value <- sum(zip_obs <= res$zip)/n_simu
  
  g1 <-ggplot(res, aes(x = zip)) +
    geom_histogram(fill = "gray", color = "black") +
    geom_vline(xintercept = zip_obs, color = "red") +
    theme_classic(base_size = 20) +
    labs(x = "Zero inflation index (expected)",
         title = "Test presence of excess zero values",
         subtitle = paste("P.value = ", p_value)) 
  print(g1)
  res
  
}
