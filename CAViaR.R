library(RcppArmadillo)
library(Rcpp)
library(numDeriv)
sourceCpp('SAVLoop.cpp')
sourceCpp('AsymmetricLoop.cpp')
sourceCpp('AdaptiveLoop.cpp')
sourceCpp('InGARCHLoop.cpp')
sourceCpp('QRObj.cpp')


CAViaR <- function(y,
                   model,
                   alpha,
                   G = 5,
                   control = list(max_iter_out = 500,
                                  max_iter_in = 500,
                                  rel_tol_QR = 1e-10,
                                  rel_tol_param = 1e-10,
                                  trace = 0),
                   lower_brent = -2,
                   upper_brent = 2,
                   seed = seed){
  
  # Correct model identification
  models <- c('SAV', 'Asymmetric', 'Adaptive', 'InGARCH')
  id_model <- which(models %in% model)
  
  # Set CAViaR_loop and CAViaR_predict as the consistent ones 
  CAViaR_loop <- get(paste(c('CAViaR', model), collapse = '_'))
  CAViaR_predict <- get(paste(c('predict_CAViaR', model), collapse = '_'))
  
  # Optimising routine constants 
  max_iter_out <- control$max_iter_out
  max_iter_in <- control$max_iter_in
  rel_tol_QR <- control$rel_tol_QR
  rel_tol_param <- control$rel_tol_param
  trace <- control$trace
  
  # Starting proposed values
  n_init <- c(10^4, 10^5, 10^4, 10^5) # As proposed by The Authors
  n <- n_init[id_model]
  
  # Consistent dimension parameter space  
  n_param <- c(3, 4, 1, 5)
  p <- n_param[id_model]
  
  # Starting points optimisation 
  n_top_val <- c(10, 15, 5, 15)
  m <- n_top_val[id_model]
  
  # Initial values
  set.seed(seed)
  init_values <- matrix(runif(p*n, 0, 0.1), ncol = p)
  
  # Empirical quantile for initializing AR-quantile sequence
  empirical_quantile <- quantile(y, alpha)
  
  if(id_model %in% c(1, 2) ){ # SAV, Asymmetric
    
    # Inital QR objective functions
    QRObj_values <- apply(init_values, 1, function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                                   data = y,
                                                                                   empiricalQuantile = empirical_quantile),
                                                            data = y,
                                                            alpha = alpha))
    # Top m beta values
    top_init_val <- init_values[order(QRObj_values, decreasing = FALSE)[1:m], ]
    top_init_QRobj <- sort(QRObj_values, decreasing = FALSE)[1:m]
    QRObj_opt_prev <- top_init_QRobj 
    
    # Optimisation 
    convergence <- 1
    j <- 1
    
    # Initializing relevant quantities
    opt_par_simpl <- matrix(NA, nrow = m, ncol = p) # Min values simplex method 
    opt_par_qn <- matrix(NA, nrow = m, ncol = p) # Min values quasi-Newton method
    
    # Optimisation loop
    while ( (convergence == 1) & (j <= max_iter_out) ) {
      
      # Simplex Algorithm
      for (i in 1:NROW(top_init_val)){
        opt <- optim(par = top_init_val[i, ],
                     fn = function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                   data = y,
                                                                   empiricalQuantile = empirical_quantile),
                                            data = y,
                                            alpha = alpha),
                     method = 'Nelder-Mead', # Simplex algorithm as fminsearch
                     control = list(trace = trace, 
                                    maxit = max_iter_in,
                                    reltol = 1e-10))
        opt_par_simpl[i, ] <- opt$par
      }
      
      # opt_par_qn <- opt_par_simpl
      
      # Quasi-Newton algorithm
      for (i in 1:NROW(opt_par_simpl)){
        opt <- optim(par = opt_par_simpl[i, ],
                     fn = function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                   data = y,
                                                                   empiricalQuantile = empirical_quantile),
                                            data = y,
                                            alpha = alpha),
                     method = 'BFGS', # Quasi-Newton as fminunc
                     control = list(trace = trace,
                                    maxit = max_iter_in,
                                    reltol = 1e-10))
        opt_par_qn[i, ] <- opt$par
      }
      
      # Convergence check
      QRObj_opt_curr <- apply(opt_par_qn, 1, function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                                      data = y,
                                                                                      empiricalQuantile = empirical_quantile),
                                                               data = y,
                                                               alpha = alpha)) # QRobj at current min values
      id_best_curr <- which.min(QRObj_opt_curr)
      
      # Relative error best vector at j-th iteration
      rel_err_param <- sum(abs(opt_par_qn[id_best_curr, ] - top_init_val[id_best_curr, ]) )/sum(abs(top_init_val[id_best_curr, ]) )
      rel_err_QRobj <- abs(QRObj_opt_curr[id_best_curr] - QRObj_opt_prev[id_best_curr])/QRObj_opt_prev[id_best_curr]
      
      if( (rel_err_param < rel_tol_param) & (rel_err_QRobj < rel_tol_QR) ){
        convergence <- 0
      } else{
        j <- j + 1 
        top_init_val <- opt_par_qn
        QRObj_opt_prev <- QRObj_opt_curr
      }
    } # End while
    
    # Output 
    opt_par_id <- which.min(QRObj_opt_curr)
    minQR <- min(QRObj_opt_curr)
    opt_beta <- opt_par_qn[opt_par_id, ]
    quantile_est <- CAViaR_loop(beta = opt_beta,
                                data = y,
                                empiricalQuantile = empirical_quantile)
    
    # Prediction via plug-in estimate
    prediction <- tail(CAViaR_predict(beta = opt_beta,
                                      data = y[length(y)],
                                      lastquantile = quantile_est[length(y)],
                                      h = 1), 1)
    
  }else if (id_model == 3){ # Adaptive
    
    # Inital QR objective functions
    QRObj_values <- apply(init_values, 1, function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                                   data = y,
                                                                                   empiricalQuantile = empirical_quantile,
                                                                                   alpha = alpha,
                                                                                   G = G),
                                                            data = y,
                                                            alpha = alpha))
    # Top m beta values
    top_init_val <- matrix(init_values[order(QRObj_values, decreasing = FALSE)[1:m], ], ncol = p)
    top_init_QRobj <- sort(QRObj_values, decreasing = FALSE)[1:m]
    QRObj_opt_prev <- top_init_QRobj 
    
    # Optimisation 
    convergence <- 1
    j <- 1
    
    # Initializing relevant quantities
    opt_par_simpl <- matrix(NA, nrow = m, ncol = p) # Min values Brent method
    opt_par_qn <- matrix(NA, nrow = m, ncol = p) # Min values quasi-Newton method
    
    # Optimisation loop
    while ( (convergence == 1) & (j <= max_iter_out) ) {
      
      # Simplex Algorithm
      for (i in 1:NROW(top_init_val)){
        opt <- optim(par = top_init_val[i, ],
                     fn = function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                   data = y,
                                                                   empiricalQuantile = empirical_quantile,
                                                                   alpha = alpha,
                                                                   G = G),
                                            data = y,
                                            alpha = alpha),
                     method = 'Brent', # No simplex because p = 1 (unidimensional optimisation problem)
                     lower = lower_brent, 
                     upper = upper_brent,
                     control = list(trace = trace, 
                                    maxit = max_iter_in,
                                    reltol = 1e-10))
        opt_par_simpl[i, ] <- opt$par
      }
      
      # Quasi-Newton algorithm
      for (i in 1:NROW(opt_par_simpl)){
        opt <- optim(par = opt_par_simpl[i, ],
                     fn = function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                   data = y,
                                                                   empiricalQuantile = empirical_quantile,
                                                                   alpha = alpha,
                                                                   G = G),
                                            data = y,
                                            alpha = alpha),
                     method = 'BFGS', # Quasi-Newton as fminunc
                     control = list(trace = trace, 
                                    maxit = max_iter_in,
                                    reltol = 1e-10))
        opt_par_qn[i, ] <- opt$par
      }
      
      # Convergence check
      QRObj_opt_curr <- apply(opt_par_qn, 1, function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                                      data = y,
                                                                                      empiricalQuantile = empirical_quantile,
                                                                                      alpha = alpha,
                                                                                      G = G),
                                                               data = y,
                                                               alpha = alpha)) # QRobj at current min values
      id_best_curr <- which.min(QRObj_opt_curr)
      
      # Relative error best vector at j-th iteration
      rel_err_param <- abs(opt_par_qn[id_best_curr] - top_init_val[id_best_curr])/abs(top_init_val[id_best_curr])
      rel_err_QRobj <- abs(QRObj_opt_curr[id_best_curr] - QRObj_opt_prev[id_best_curr])/QRObj_opt_prev[id_best_curr]
      
      if((rel_err_param < rel_tol_param) & (rel_err_QRobj < rel_tol_QR)){
        convergence <- 0
      } else{
        j <- j + 1 
        top_init_val <- opt_par_qn
        QRObj_opt_prev <- QRObj_opt_curr
      }
    }
    
    # Output
    minQR <- min(QRObj_opt_curr)
    opt_par_id <- which.min(QRObj_opt_curr)
    opt_beta <- opt_par_qn[opt_par_id, ]
    quantile_est <- CAViaR_loop(beta = opt_beta,
                                data = y,
                                empiricalQuantile = empirical_quantile,
                                alpha = alpha,
                                G = G)
    
    # Prediction
    prediction <- tail(CAViaR_predict(beta = opt_beta,
                                      data = y[length(y)],
                                      lastquantile = quantile_est[length(y)],
                                      alpha = alpha,
                                      G = G,
                                      h = 1), 1)
  } else{ # InGARCH
    
    # Inital QR objective functions
    QRObj_values <- apply(init_values, 1, function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                                   data = y,
                                                                                   empiricalQuantile = empirical_quantile,
                                                                                   alpha = alpha),
                                                            data = y,
                                                            alpha = alpha))
    # Top m beta values
    top_init_val <- init_values[order(QRObj_values, decreasing = FALSE)[1:m], ]
    top_init_QRobj <- sort(QRObj_values, decreasing = FALSE)[1:m]
    QRObj_opt_prev <- top_init_QRobj 
    
    # Optimisation 
    convergence <- 1
    j <- 1
    
    # Initializing relevant quantities
    opt_par_simpl <- matrix(NA, nrow = m, ncol = p) # Min values simplex method 
    opt_par_qn <- matrix(NA, nrow = m, ncol = p) # Min values quasi-Newton method
    
    # Optimisation loop
    while ( (convergence == 1) & (j <= max_iter_out) ) {
      
      # Simplex Algorithm
      for (i in 1:NROW(top_init_val)){
        opt <- optim(par = top_init_val[i, ],
                     fn = function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                   data = y,
                                                                   empiricalQuantile = empirical_quantile,
                                                                   alpha = alpha),
                                            data = y,
                                            alpha = alpha),
                     method = 'Nelder-Mead', # Simplex algorithm as fminsearch
                     control = list(trace = trace,
                                    maxit = max_iter_in,
                                    reltol = 1e-10))
        opt_par_simpl[i, ] <- opt$par
      }
      
      # Quasi-Newton algorithm
      for (i in 1:NROW(opt_par_simpl)){
        opt <- optim(par = opt_par_simpl[i, ],
                     fn = function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                   data = y,
                                                                   empiricalQuantile = empirical_quantile,
                                                                   alpha = alpha),
                                            data = y,
                                            alpha = alpha),
                     method = 'L-BFGS-B', # Quasi-Newton as fminunc
                     lower = c(-Inf, rep(1e-3, 4)),
                     upper = c(rep(Inf, 5)),
                     control = list(trace = trace,
                                    maxit = max_iter_in),
                     gr = gr)
        print(opt$convergence)
        opt_par_qn[i, ] <- opt$par
      }
      
      # Convergence check
      QRObj_opt_curr <- apply(opt_par_qn, 1, function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                                                      data = y,
                                                                                      empiricalQuantile = empirical_quantile,
                                                                                      alpha = alpha),
                                                               data = y,
                                                               alpha = alpha)) # QRobj at current min values
      id_best_curr <- which.min(QRObj_opt_curr)
      
      # Relative error best vector at j-th iteration
      rel_err_param <- sum(abs(opt_par_qn[id_best_curr, ] - top_init_val[id_best_curr, ]) )/sum(abs(top_init_val[id_best_curr, ]) )
      rel_err_QRobj <- abs(QRObj_opt_curr[id_best_curr] - QRObj_opt_prev[id_best_curr])/QRObj_opt_prev[id_best_curr]
      
      if( (rel_err_param < rel_tol_param) & (rel_err_QRobj < rel_tol_QR) ){
        convergence <- 0
      } else{
        j <- j + 1 
        top_init_val <- opt_par_qn
        QRObj_opt_prev <- QRObj_opt_curr
      }
    } # End while
    
    # Output 
    opt_par_id <- which.min(QRObj_opt_curr)
    minQR <- min(QRObj_opt_curr)
    opt_beta <- opt_par_qn[opt_par_id, ]
    quantile_est <- CAViaR_loop(beta = opt_beta,
                                data = y,
                                empiricalQuantile = empirical_quantile,
                                alpha = alpha)
    
    # Prediction via plug-in estimate
    prediction <- tail(CAViaR_predict(beta = opt_beta,
                                      data = y[length(y)],
                                      lastquantile = quantile_est[length(y)],
                                      h = 1), 1)
  }
  
  # Output
  ret <- list(alpha = alpha,
              quantile = quantile_est,
              opt_beta = opt_beta,
              minQR = minQR,
              prediction = prediction,
              convergence = convergence)
  
  return(ret)
}







# gr <- function(par){
#   numDeriv::grad(function(x) QRObj(Quantile = CAViaR_loop(beta = x,
#                                                           data = y,
#                                                           empiricalQuantile = empirical_quantile,
#                                                           alpha = alpha),
#                                    data = y,
#                                    alpha = alpha),
#                  x = par)
# }