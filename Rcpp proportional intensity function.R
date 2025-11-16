library(Rcpp)

# ---------------------------------------------------------
# RCPP IMPLEMENTATIONS
# ---------------------------------------------------------
Rcpp::cppFunction('
double log_likelihood_cpp(const NumericMatrix& X,
                          const NumericMatrix& N_ij,
                          const NumericMatrix& I_ij,
                          const NumericVector& beta,
                          const NumericVector& Lambda0,
                          const NumericVector& W) {
  int N = X.nrow();
  int p = X.ncol();
  int J = Lambda0.size();
  double ll = 0.0;

  for (int i = 0; i < N; ++i) {
    double xb = 0.0;
    for (int k = 0; k < p; ++k)
      xb += beta[k] * X(i, k);
    double exp_xb = std::exp(xb);

    for (int j = 0; j < J; ++j) {
      if (I_ij(i, j) == 0) continue;
      double mu = W[i] * Lambda0[j] * exp_xb;
      double y = N_ij(i, j);
      if (mu > 0.0) {
        ll += y * std::log(mu) - mu - std::lgamma(y + 1.0);
      }
    }
  }
  return ll;
}
')

Rcpp::cppFunction('
NumericVector update_Lambda0_cpp(const NumericMatrix& X,
                                 const NumericMatrix& N_ij,
                                 const NumericMatrix& I_ij,
                                 const NumericVector& beta,
                                 const NumericVector& W,
                                 double shape_prior = 0.05,
                                 double rate_prior = 1.0) {
  int N = X.nrow();
  int p = X.ncol();
  int J = I_ij.ncol();
  NumericVector Lambda0(J);

  for (int j = 0; j < J; ++j) {
    double shape_post = shape_prior;
    double rate_post = rate_prior;

    for (int i = 0; i < N; ++i) {
      if (I_ij(i, j) == 1) {
        double xb = 0.0;
        for (int k = 0; k < p; ++k) xb += beta[k] * X(i, k);
        double exp_xb = std::exp(xb);
        shape_post += N_ij(i, j);
        rate_post += W[i] * exp_xb;
      }
    }
    Lambda0[j] = R::rgamma(shape_post, 1.0 / rate_post);
  }

  return Lambda0;
}
')

Rcpp::cppFunction('
NumericVector update_W_cpp(const NumericMatrix& X,
                           const NumericMatrix& N_ij,
                           const NumericMatrix& I_ij,
                           const NumericVector& beta,
                           const NumericVector& Lambda0,
                           double eta) {
  int N = X.nrow();
  int p = X.ncol();
  int J = Lambda0.size();
  NumericVector W(N);

  for (int i = 0; i < N; ++i) {
    double xb = 0.0;
    for (int k = 0; k < p; ++k) xb += beta[k] * X(i, k);
    double exp_xb = std::exp(xb);

    double count_sum = 0.0;
    double exposure_sum = 0.0;
    for (int j = 0; j < J; ++j) {
      if (I_ij(i, j) == 1) {
        count_sum += N_ij(i, j);
        exposure_sum += Lambda0[j] * exp_xb;
      }
    }

    double shape_post = eta + count_sum;
    double rate_post = eta + exposure_sum;
    W[i] = R::rgamma(shape_post, 1.0 / rate_post);
  }

  return W;
}
')

# ---------------------------------------------------------
# MAIN FUNCTION (R)
# ---------------------------------------------------------
fit_proportional_intensity <- function(X_train,
                                       X_test,
                                       recurrent_train,
                                       recurrent_test,
                                       terminal_train,
                                       terminal_test,
                                       num_thin = 1,
                                       num_burn,
                                       num_save) {
  # --- Construct time grid ---
  all_list <- c(recurrent_train, recurrent_test)
  all_recurrence_time <- unique(sort(unlist(all_list)))
  J <- length(all_recurrence_time)
  
  N_train <- nrow(X_train)
  p <- ncol(X_train)
  
  # --- Initialize N_ij and I_ij ---
  N_ij <- array(0, dim = c(N_train, J))
  I_ij <- array(0, dim = c(N_train, J))
  
  for (i in 1:N_train) {
    recurrent <- recurrent_train[[i]]
    if (length(recurrent) == 0) next
    k <- match(recurrent, all_recurrence_time)
    k <- k[!is.na(k)]
    if (length(k) > 0) N_ij[i, k] <- 1
  }
  for (i in 1:N_train) {
    terminal <- terminal_train[i]
    if (is.na(terminal) || length(terminal) == 0) next
    k <- which(all_recurrence_time <= terminal)
    if (length(k) > 0) I_ij[i, k] <- 1
  }
  
  # --- MCMC storage ---
  num_iter <- num_burn + num_thin * num_save
  beta_samples    <- array(NA, dim = c(num_save, p))
  Lambda0_samples <- array(NA, dim = c(num_save, J))
  W_samples       <- array(NA, dim = c(num_save, N_train))
  eta_samples     <- numeric(num_save)
  
  # --- Initial values ---
  beta <- rep(0, p)
  Lambda0 <- rep(1, J)
  W <- rep(1, N_train)
  eta <- 40
  
  log_prior_beta <- function(b) sum(dnorm(b, mean = 0, sd = 1, log = TRUE))
  
  iter_save <- 0
  for (iter in 1:num_iter) {
    
    # --- Update beta (Random-Walk MH) ---
    beta_prop <- beta + rnorm(p, 0, 0.1)
    log_post_curr <- log_likelihood_cpp(X_train, N_ij, I_ij, beta, Lambda0, W) + log_prior_beta(beta)
    log_post_prop <- log_likelihood_cpp(X_train, N_ij, I_ij, beta_prop, Lambda0, W) + log_prior_beta(beta_prop)
    acc <- log_post_prop - log_post_curr
    if (is.finite(acc) && log(runif(1)) < acc) beta <- beta_prop
    
    # --- Update Lambda0 and W (Rcpp Gibbs steps) ---
    Lambda0 <- update_Lambda0_cpp(X_train, N_ij, I_ij, beta, W)
    W <- update_W_cpp(X_train, N_ij, I_ij, beta, Lambda0, eta)
    
    # --- Update eta using diversitree::mcmc ---
    loglik_eta_given_W <- function(eta_val) {
      if (eta_val <= 0) return(-Inf)
      sum(dgamma(W, shape = eta_val, rate = eta_val, log = TRUE))
    }
    
    eta_update <- diversitree::mcmc(
      lik    = loglik_eta_given_W,
      nsteps = 1,
      w      = 1,
      x.init = eta,
      prior  = function(x) dgamma(x, shape = 40, rate = 0.1, log = TRUE),
      lower  = 20,
      upper  = Inf
    )
    eta <- eta_update$pars
    
    # --- Save samples ---
    if (iter > num_burn && (iter - num_burn) %% num_thin == 0) {
      iter_save <- iter_save + 1
      beta_samples[iter_save, ]    <- beta
      Lambda0_samples[iter_save, ] <- Lambda0
      W_samples[iter_save, ]       <- W
      eta_samples[iter_save]       <- eta
    }
    
    # --- Progress ---
    if (iter %% 100 == 0) cat("Iter:", iter, "\n")
  }
  
  # --- Return samples ---
  list(
    beta = beta_samples,
    Lambda0 = Lambda0_samples,
    W = W_samples,
    eta = eta_samples
  )
}
