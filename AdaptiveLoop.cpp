#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector CAViaR_Adaptive(NumericVector beta, NumericVector data, double alpha, double empiricalQuantile, double G)
{
  int i;
  int data_length = data.length();
  NumericVector Quantile(data_length);
  
  // Initialize output variables
  Quantile[0] = empiricalQuantile;
  
  for(i = 1; i < data_length; i++)
  {
    // Adaptive 
    Quantile[i] = Quantile[i - 1] + beta[0]*( pow(1 + exp(G*(data[i - 1] - Quantile[i - 1]) ), -1.0) - alpha);
  }
  return(Quantile);
}



// [[Rcpp::export]]
NumericVector predict_CAViaR_Adaptive(NumericVector beta, NumericVector data, double alpha, double lastquantile, double G, int h)
{
  int i;
  int data_length = data.length();
  NumericVector Quantile(data_length + h);
  
  // Initialize output variables
  Quantile[0] = lastquantile;
  
  for(i = 1; i < h + 1; i++)
  {
    // Adaptive 
    Quantile[i] = Quantile[i - 1] + beta[0]*( pow(1 + exp(G*(data[i - 1] - Quantile[i - 1]) ), -1.0) - alpha);
  }
  return(Quantile);
}