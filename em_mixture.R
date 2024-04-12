library(MASS)
library(mvtnorm)
library(tidyr)
library(ggplot2)

# unsupervised clustering algorithm via EM
# we'll start with a bivariate, 2-component mixture of normals and then extend to more components and dimensions

# -------------E step-------------
# priors, mus, sigmas refer to parameter choice in current iteration
calulate_posterior <- function(i, j, data, priors, mus, sigmas, g) {
  # tupel of data (p-dimensional vector)
  data_point <- data[i,]
  prior <- priors[j]
  mu <- mus[j,]
  sigma <- sigmas[[j]]
  denominator <- 0
  # hier wird component durch seq_len(anzahlAnKlassen) zählen
  # vlt kann man das noch effizienter machen anstatt mit for-Loop
  for (component in seq_len(g)) {
    denominator <- denominator + (priors[component] * (exp(-0.5 * t(data_point - mus[component,]) %*% solve(sigmas[[component]]) %*% (data_point - mus[component,])) / sqrt((2 * pi)^2 * det(sigmas[[component]]))))
  }
  return((prior * (exp(-0.5 * t(data_point - mu) %*% solve(sigma) %*% (data_point - mu)) / sqrt((2 * pi)^2 * det(sigma)))) / denominator)
}
# --------------------------------

# -------------M step-------------
# estimate priors
estimate_prior <- function(data, j, priors, mus, sigma, g) {
  priors_j_old <- priors
  mu_j_old <- mus
  sigma_j_old <- sigma
  post_record <- vector()
  # for durch apply ersetzen, dann returnt die Funktion aber keinen Vektor sondern schon die Summe, deswegen muss das Produkt dann angepasst werden
  for (i in seq_len(nrow(data))) {
    post_j <- calulate_posterior(i, j, data, priors_j_old, mu_j_old, sigma_j_old, g)
    post_record <- c(post_record, post_j)
  }
  return(post_record)
}

# estimate mus
estimate_mu <- function(data, j, sum_post, priors, mus, sigma, g) {
  priors_j_old <- priors
  mu_j_old <- mus
  sigma_j_old <- sigma
  weights_record <- vector()
  sum_post <- sum_post
  for (i in seq_len(nrow(data))) {
    # hier statt for Schleife apply(Produkt und dann das aufsummieren)
    weights_j <- (calulate_posterior(i, j, data, priors_j_old, mu_j_old, sigma_j_old, g) / sum_post) * data[i,]
    weights_record <- rbind(weights_record, weights_j)
  }
  return(apply(weights_record, MARGIN = 2, sum))
}

# estimate sigmas
estimate_sigma <- function(data, j, sum_post, mu_j_new, priors, mus, sigma, g) {
  components_record <- list()
  sum_post <- sum_post
  mu_j_new <- mu_j_new
  priors_j_old <- priors
  mu_j_old <- mus
  sigma_j_old <- sigma
  for (i in seq_len(nrow(data))) {
    weight <- (calulate_posterior(i, j, data, priors_j_old, mu_j_old, sigma_j_old, g) / sum_post)
    # in der Matrix muss nrow und ncol dann variabel sein
    component_j <- diag(as.numeric(weight), 2) %*% ((data[i,] - mu_j_new) %*% t(data[i,] - mu_j_new))
    components_record[[i]] <- component_j
  }
  return(Reduce("+", components_record))
}
# --------------------------------

# generate data from true parameters
# ------- 1 (2 components)
n <- 1000
mu1 <- c(-1, -1)
mu2 <- c(1, 1)
mus <- rbind(mu1,mu2)
sigma1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
sigma2 <- matrix(c(0.24, 0.5, 0.5, 8), nrow = 2)
sigmas <- list(sigma1,sigma2)
n1 <- 600
n2 <- 400
x1 <- rmvnorm(n1, mu1, sigma1)
x2 <- rmvnorm(n2, mu2, sigma2)
data <- rbind(x1, x2)

# ------- 2 (3 components)
n <- 1000
mu1 <- c(-10, -10)
mu2 <- c(3, 3)
mu3 <- c(10, 10)
mus <- rbind(mu1,mu2,mu3)
sigma1 <- matrix(c(1, 0, 0, 1), nrow = 2)
sigma2 <- matrix(c(1, 0, 0, 1), nrow = 2)
sigma3 <- matrix(c(1, 0, 0, 1), nrow = 2)
sigmas <- list(sigma1,sigma2,sigma3)
n1 <- 333
n2 <- 333
n3 <- 334
x1 <- rmvnorm(n1, mu1, sigma1)
x2 <- rmvnorm(n2, mu2, sigma2)
x3 <- rmvnorm(n3, mu3, sigma3)
data <- rbind(x1, x2, x3)

# create data set with higher overlap and three components
# ------- 3 (3 components)
n <- 1000
mu1 <- c(-3, -3)
mu2 <- c(1, 1)
mu3 <- c(2, -1)
mus <- rbind(mu1,mu2,mu3)
sigma1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
sigma2 <- matrix(c(0.24, 0.5, 0.5, 8), nrow = 2)
sigma3 <- matrix(c(2, 1, 1, 2), nrow = 2)
sigmas <- list(sigma1,sigma2,sigma3)
n1 <- 100
n2 <- 450
n3 <- 350
x1 <- rmvnorm(n1, mu1, sigma1)
x2 <- rmvnorm(n2, mu2, sigma2)
x3 <- rmvnorm(n3, mu3, sigma3)
data <- rbind(x1, x2, x3)

# ------- init 3 (3 components)
# initiate random starting points
mu1 <- c(0.1, 0)
mu2 <- c(0.2, 0)
mu3 <- c(0.3, 0)
mu_j_old <- rbind(mu1,mu2,mu3)
sigma1 <- matrix(c(1, 0, 0, 1), nrow = 2)
sigma2 <- matrix(c(1, 0, 0, 1), nrow = 2)
sigma3 <- matrix(c(1, 0, 0, 1), nrow = 2)
sigma_j_old <- list(sigma1,sigma2,sigma3)
priors_j_old <- c(0.05, 0.45, 0.5)

em <- function(data, priors_j_old, mu_j_old, sigma_j_old, g) {
  previous_mus <- 0
  for (i in seq_len(100)) {
    # data grid has to be made variable too later but we'll keep 3 for now
    data.grid <- expand.grid(s.1 = seq(-50, 50, length.out=200), 
                             s.2 = seq(-50, 50, length.out=200))
    
    geom_list <- list()
    
    # Erstelle Geometrieobjekte für jedes j und füge sie zur Liste hinzu
    for (j in seq_len(g)) {
      q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu_j_old[j,], sigma = sigma_j_old[[j]]))
      geom_list[[j]] <- geom_contour(data = q.samp, aes(x = s.1, y = s.2, z = prob), color = "darkred", alpha = 0.5)
    }
    
    # Erstelle ggplot mit den Punkten und den Konturlinien für jedes j
    print(ggplot(as.data.frame(data), aes(x = V1, y = V2)) +
      geom_point(alpha = 0.5) +
      geom_list +
      theme_minimal())
    
    
    priors_vector <- vector()
    mus_matrix <- matrix(data = NA, nrow = g, ncol = 2)
    sigma_list <- list()
    for (j in seq_len(g)) {
      post_rec <- estimate_prior(data, j, priors_j_old, mu_j_old, sigma_j_old, g)
      sum_post <- sum(post_rec)
      prior_j_new <- mean(post_rec)
      priors_j_new <- mean(estimate_prior(data, j, priors_j_old, mu_j_old, sigma_j_old, g))
      mu_j_new <- estimate_mu(data, j, sum_post, priors_j_old, mu_j_old, sigma_j_old, g)
      sigma_j_new <- estimate_sigma(data, j, sum_post, mu_j_new, priors_j_old, mu_j_old, sigma_j_old, g)
      
      priors_vector <- c(priors_vector, priors_j_new)
      mus_matrix[j,] <- mu_j_new
      sigma_list[[j]] <- sigma_j_new
    }
    previous_mus <- mu_j_old[1,1]
    priors_j_old <- priors_vector
    mu_j_old <- mus_matrix
    sigma_j_old <- sigma_list
    print(mus_matrix)
    
  }
}

# run the following line to try the em function in an isolated way (dataset is created above)
# em(data_matrix, priors_j_old, mu_j_old, sigma_j_old, 3)

# label data for illustration with ggplot
data_with_labels <- data.frame(data)
label_vector <- rep(1, n1)
label_vector <- c(label_vector, rep(2, n2))
label_vector <- c(label_vector, rep(3, n3))
data_with_labels$label <- label_vector

# this function generates data with sample size n, some number of classes (hidden
# components), a certain bayes error (overlap, not yet implemented yet), and a 
# certain correlation between the variables. The prior ratio can only be set to "equal"
# meaning all classes have similar prior probabilities (each sampled from a normal
# distribution centered around 1 / number_of_classes)
generate_data <- function(n, number_of_classes, bayes_error, correlation, prior_ratio = "equal") {
  mu_matrix <- matrix(data = NA, nrow = number_of_classes, ncol = 2)
  sigma_list <- list()
  prior_vector <- vector(mode = "numeric", length = number_of_classes)
  max <- 1
  x <- matrix(data = NA, nrow = 0, ncol = 3)
  
  for (i in seq_len(number_of_classes - 1)) {
    mu_matrix[i, ] <- runif(2, 0, 10)
    
    current_sigma <- matrix(c(1, 0, 0, 1), nrow = 2)
    sigma_list[[i]] <- current_sigma
    
    ratio_accepted <- FALSE
    while(ratio_accepted == FALSE) {
      ratio <- rnorm(1, 1 / number_of_classes, 0.01)
      if (ratio > max) {
        ratio_accepted <- FALSE
      } else {
        ratio_accepted <- TRUE
      }
    }
    prior_vector[i] <- ratio
    max <- max - ratio
    
    x_new <- rmvnorm(round(n * ratio), mu_matrix[i, ], current_sigma)
    labels <- rep(i, round(n * ratio))
    x_new <- cbind(x_new, labels)
    x <- rbind(x, x_new)
    x_new <- x
    
    
  }
  
  mu_matrix[i + 1, ] <- runif(2, 0, 10)
  current_sigma <- matrix(c(1, 0, 0, 1), nrow = 2)
  sigma_list[[i + 1]] <- current_sigma
  
  prior_vector[length(prior_vector)] <- max
  
  
  x_new <- rmvnorm(round(n * max), mu_matrix[i + 1, ], current_sigma)
  labels <- rep(i + 1, round(n * max))
  x_new <- cbind(x_new, labels)
  x <- rbind(x, x_new)
  x_new <- x
  
  return(x)
}

# this now combines random data generation with the EM algorithm
run_algorithm <- function(sample_size, amount_of_distributions = 2) {
  data_df <- as.data.frame(generate_data(sample_size, amount_of_distributions, 1, 1, "equal"))
  print(ggplot(data = as.data.frame(data_df), aes(x = V1, y = V2, color = as.factor(labels))) +
    geom_point() +
    theme_minimal())
  # wait for 3 sec
  Sys.sleep(3)
  data_matrix <- as.matrix(data_df[,-3])
  em(data_matrix, priors_j_old, mu_j_old, sigma_j_old, amount_of_distributions)
  ggplot(data = as.data.frame(data_df), aes(x = V1, y = V2, color = as.factor(labels))) +
    geom_point() +
    theme_minimal()
}

# test run
run_algorithm(2000, 3)

# also try with iris data
data <- iris
# only take "sepal.length", "petal.width"
data <- data[,c(2,4,5)]
ggplot(data = data, aes(x = Sepal.Width, y = Petal.Width, color = as.factor(Species))) +
  geom_point() +
  theme_minimal()

iris_matrix <- as.matrix(data[,-3])
colnames(iris_matrix) <- c("V1", "V2")

# run em
em(iris_matrix, priors_j_old, mu_j_old, sigma_j_old, 3)



set.seed(12)
set.seed(13)
# doesn't work with arg 1 yet
data_df <- as.data.frame(generate_data(2000, 3, 1, 1, "equal"))
data_matrix <- as.matrix(data_df[,-3])


# mit Farben
ggplot(data = as.data.frame(data_df), aes(x = V1, y = V2, color = as.factor(labels))) +
  geom_point() +
  theme_minimal()

# ohne Farben
ggplot(data = data_df, aes(x = X1, y = X2)) +
  geom_point() +
  theme_minimal()
  









