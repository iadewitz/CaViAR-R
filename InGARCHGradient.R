grad <- function(par){
  numDeriv::grad(function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                          data = y,
                                                          empiricalQuantile = empirical_quantile),
                                   data = y,
                                   alpha = alpha),
                 x = par)
}







numDeriv::grad(function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                        data = y,
                                                        empiricalQuantile = empirical_quantile,
                                                        alpha = alpha),
                                 data = y,
                                 alpha = alpha),
               x = opt_beta)


opt1 <- optim(par = opt_par_simpl[i, ],
             fn = function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                           data = y,
                                                           empiricalQuantile = empirical_quantile),
                                    data = y,
                                    alpha = alpha),
             method = 'BFGS', # Quasi-Newton as fminunc
             gr = grad,
             control = list(trace = trace,
                            maxit = max_iter_in,
                            reltol = 1e-10))










numDeriv::hessian(function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                           data = y,
                                                           empiricalQuantile = empirical_quantile,
                                                           alpha = alpha),
                                    data = y,
                                    alpha = alpha),
                  x = abs(opt_beta))




opt_beta %*% (numDeriv::hessian(function(x) QRObj(Quantile = CAViaR_loop(beta = x,
                                                           data = y,
                                                           empiricalQuantile = empirical_quantile),
                                    data = y,
                                    alpha = alpha),
                  x = opt_beta) ) %*% opt_beta

