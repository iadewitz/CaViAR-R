#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sgn(double v) {
  double sign;
  sign = 1*(v >= 0) + -1*(v < 0);
  return(sign);
}

// [[Rcpp::export]]
NumericVector CAViaR_InGARCH(NumericVector beta, NumericVector data, double empiricalQuantile, double alpha) {
  int i;
  int data_length = data.length();
  NumericVector Quantile(data_length);
  
  // Initialize output variables 
  Quantile[0] = empiricalQuantile;
  
  for(i = 1; i < data_length; i++)  {
    // Indirect GARCH(1, 1)
    Quantile[i] = beta[0] + sgn(alpha - 0.5)*sqrt(beta[1]*(beta[2] + beta[3]*pow(data[i - 1] - beta[0], 2) ) + beta[4]*pow(Quantile[i - 1] - beta[0], 2) );
  }
  return(Quantile);
}



// [[Rcpp::export]]
NumericVector predict_CAViaR_InGARCH(NumericVector beta, NumericVector data, double lastquantile, int h){
  int i;
  int data_length = data.length();
  NumericVector Quantile(data_length + h);
  
  // Initialize output variables
  Quantile[0] = lastquantile;
  
  for(i = 1; i < h + 1; i++){
    // Indirect GARCH(1, 1)
    Quantile[i] = beta[0] + sgn(lastquantile)*sqrt(beta[1]*(beta[2] + beta[3]*pow(data[i - 1] - beta[0], 2) ) + beta[4]*pow(Quantile[i - 1] - beta[0], 2) );
  }
  return(Quantile);
}


