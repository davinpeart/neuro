negative_binomial <- function(predict_ancillary, repeated_measures, priors) {
  
  # generate stancode ==========================================================
  # functions block
  functions <- "
  "
  
  # data block
  data <- paste0(
    "
// observations and fixed effects
    int<lower=1> N;  // number of observations
    int Y[N];  // response variable
    int<lower=1> K;  // number of fixed effects
    matrix[N, K] X;  // fixed effect design matrix", 
    if(predict_ancillary) {"
// effects on distributional parameter(s)
    int<lower=1> K_phi;  // number of effects in phi design
    matrix[N, K_phi] X_phi;  // phi design matrix"},
    if(repeated_measures) {
  "
// random effects
    int<lower=1> N_I;  // number of subjects
    int<lower=1,upper=N_I> I[N];  // subject identifier
    int<lower=1> J;  // number of distributional parameters
"
}
  )
  
  # parameters block
  parameters <- paste0(
    "
// fixed effects
    vector[K] beta;  // regression coefficients
    ", if(!predict_ancillary) "real phi" else "vector[K_phi] beta_phi",
    ";  // log positive phi
     ",
    if(repeated_measures) {
      "
// random effects
    matrix[J, N_I] z_I;  // standardized subject intercepts
    vector<lower=0>[J] sigma_I;  // sd for subject intercepts and slopes
    cholesky_factor_corr[J] L_I;  // correlation matrix for subject intercepts and slopes
"
    }
  )
  
  # print transformed parameters
  transformed_parameters <-
    if(repeated_measures) {
      "
// random effects
    matrix[J, N_I] z; // non-centered subject intercepts and slopes
    z = diag_pre_multiply(sigma_I, L_I) * z_I;
"
    } else { "
      "}
  
  # print priors
  priors <- paste0("
// fixed effect priors
    beta ~ ", priors[["beta"]], ";
    ", if(predict_ancillary) {
    "beta_"
  },
  "phi ~ ", priors[["phi"]], ";
  ", 
  if(repeated_measures) {
    paste0(
      "
// random effect priors
    L_I ~ ", priors[["cor"]], ";
    sigma_I ~ ", priors[["sd"]], ";
    to_vector(z_I) ~ std_normal();  // standard normal on standardized effects
"
    )
  }
  )

  # print likelihood
  mu_exp <- paste0("X[n, 1] * (beta[1]", if(repeated_measures) " + z[1, I[n]]",
                       ") + X[n, 2:K] * beta[2:K]", collapse = "")
  if(!predict_ancillary) {
    phi_exp <- paste0("exp(phi", if(repeated_measures) " + z[2, I[n]]", ")")
  } else {
    phi_exp <- paste0("exp(X_phi[n, 1] * (beta_phi[1]", if(repeated_measures) " + z[2, I[n]]",
                      ") + X_phi[n, 2:K_phi] * beta_phi[2:K_phi])", collapse = "")
  }
  
  likelihood <- paste0("
// likelihood
    for(n in 1:N) { 
      target += neg_binomial_2_log_lpmf(Y[n] | "
      , mu_exp, ", 
      ", phi_exp, ");
    }
"
  )
  
  # print generated quantities
  generated_quantities <- paste0(
    if(repeated_measures) {"
// recover omega
    matrix[J, J] omega;
    omega = multiply_lower_tri_self_transpose(L_I);

// store random intercepts as array of vectors
    array[J] vector[N_I] z_array;
    for(j in 1:J) {
      z_array[j] = to_vector(z[j, ]);
    }

// correlated replications of random intercepts
    array[J] vector[N_I] z_rep;  // random intercept replications
    z_rep = multi_normal_rng(z_array, quad_form_diag(omega, sigma_I));
      "
    }, "
// replications for posterior predictive checks",
    paste0("
    array[N] real y_rep; // Y replications
    for(n in 1:N) {
        y_rep[n] = neg_binomial_2_log_rng(
        ", mu_exp, ", 
        ", phi_exp, ");
    }
"
    )
    )
  
  stan_code <- paste0(
"functions {", functions, "}
data {", data, "}
parameters {", parameters, "}
transformed parameters {", transformed_parameters, "}
model {", priors, likelihood, "}
generated quantities {", generated_quantities, "}")
  
  return(list(functions = functions, data = data, parameters = parameters,
              transformed_parameters = transformed_parameters,
              priors = priors, likelihood = likelihood,
              generated_quantities = generated_quantities,
              stan_code = stan_code))
}

attr(negative_binomial, "npar") <- 2
attr(negative_binomial, "parnames") <- c("mu", "phi")

fit_model <-
  function(data, dv, iv, id = NULL, family, dpar_iv = NULL, priors, 
           chains = 4, warmup = 5000, iter = 10000, cores = 4, refresh = 2,
           seed = runif(n = 1, min = 1, max = 99999), control = list(
             adapt_delta = .99, stepsize = .3, max_treedepth = 15), 
           standardize_coef = FALSE, stancode_only = FALSE, init = NULL) {
    
    # model matrices with treatment contrasts ==================================
    # list of design matrices for distributional parameters
    parnames <- attr(family, "parnames")
    npar <- attr(family, "npar")
    list_mat <- vector("list", npar)
    names(list_mat) <- parnames

    # design matrix for mu
    for(i in 1:length(iv)) {
      attr(data[[iv[i]]], "contrasts") <- 
        contr.treatment(levels(data[[iv[i]]]))
    }
    list_mat[["mu"]] <- model.matrix(formula(paste0(
      dv, "~", "1+", paste0(iv, collapse = "*"))), data)

    # design matrices for all other distributional parameters
    if(!is.null(dpar_iv)) {
      for(i in 1:length(dpar_iv)) {
        for(n in 1:length(dpar_iv[[i]])) {
          attr(data[[dpar_iv[[i]][n]]], "contrasts") <- 
            contr.treatment(levels(data[[dpar_iv[[i]][n]]]))
        }
      }
      for(i in 1:length(dpar_iv)) {
        list_mat[[names(dpar_iv[i])]] <- model.matrix(formula(paste0(
          dv, "~", "1+", paste0(dpar_iv[[i]], collapse = "*"))), data)
      }
    }

    # standardize response variable if required ================================
    if(standardize_coef) {
      data[[dv]] <- data[[dv]] / sd(data[[dv]])
    }
    
    # generate stancode ========================================================
    family_list <- family(predict_ancillary = !is.null(dpar_iv),
                          repeated_measures = !is.null(id),
                          priors = priors)
    stan_code <- family_list[["stan_code"]]
    
    # return stancode without fitting model if required ========================
    if(stancode_only) {
      return(cat(stan_code))
    }
    
    # create stan data =========================================================
    stan_data <- list(
      N = length(data[[dv]]),
      Y = data[[dv]],
      K = ncol(list_mat[["mu"]]),
      X = list_mat[["mu"]]
    )
    
    if(!is.null(id)) {
      stan_data[["I"]] <- as.integer(data[[id]])
      stan_data[["N_I"]] <- length(levels(data[[id]]))
      stan_data[["J"]] <- npar
    }
    
    if(!is.null(dpar_iv)) {
      for(i in 2:npar) {
        stan_data[[paste("X_", names(list_mat[i]), sep = "")]] <- list_mat[[i]]
        stan_data[[paste("K_", names(list_mat[i]), sep = "")]] <- ncol(list_mat[[i]])
      }
    }
    
    # fit stan model ===========================================================
    model <-
      rstan::sampling(
        object = rstan::stan_model(
          model_code = stan_code),
        data = stan_data,
        chains = chains,
        warmup = warmup,
        iter = iter,
        cores = cores,
        refresh = refresh,
        init = if(!is.null(init)) init else "random",
        seed = seed,
        control = control
    )
    
    # set names ================================================================
    names(model)[1:ncol(list_mat[["mu"]])] <- colnames(list_mat[["mu"]])
    
    return(model)
  }




