#include <Rcpp.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double QRObj(NumericVector Quantile, NumericVector data, double alpha)
{
  int i;
  int data_length = data.length();
  NumericVector QRObj_ind(data_length);
  
  for(i = 0; i < data_length; i++)
  {
    // Objective function Regression Quantile
    QRObj_ind[i] = (1.0/data_length)*(alpha - (data[i] < Quantile[i]))*(data[i] - Quantile[i]);
  }
  
  double QRObj_sum = accumulate(QRObj_ind.begin(), QRObj_ind.end(), 0.0);

  return(QRObj_sum);
}
