#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector CAViaR_InEGARCH(NumericVector beta, NumericVector data, double empiricalQuantile)
{
  int i;
  int data_length = data.length();
  NumericVector Quantile(data_length);
  
  // Initialize output variables
  Quantile[0] = empiricalQuantile;
  
  for(i = 1; i < data_length; i++)
  {
    // Indirect EGARCH (1, 1)
    Quantile[i] = beta[0] + beta[1] * Quantile[i - 1] + beta[2]*(data[i - 1]*(data[i - 1] >= 0)) + beta[3]*(-data[i - 1]*(data[i - 1] < 0) );
  }
  return(Quantile);
}



// [[Rcpp::export]]
NumericVector predict_CAViaR_InEGARCH(NumericVector beta, NumericVector data, double lastquantile, int h)
{
  int i;
  int data_length = data.length();
  NumericVector Quantile(data_length + h);
  
  // Initialize output variables
  Quantile[0] = lastquantile;
  
  for(i = 1; i < h + 1; i++)
  {
    // Indirect EGARCH (1, 1)
    Quantile[i] = beta[0] + beta[1] * Quantile[i - 1] + beta[2]*(data[i - 1]*(data[i - 1] >= 0)) + beta[3]*(-data[i - 1]*(data[i - 1] < 0));
  }
  return(Quantile);
}