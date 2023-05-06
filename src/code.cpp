#include <Rcpp.h>
using namespace Rcpp;

int mod(int x, int y) {
  int r = x % y;
  return r < 0 ? r + y : r;
}

// [[Rcpp::export]]
double energy_diff3(double focal, double neigh, double B, double J, double inertia) {
  double en_diff = 2 * (B + neigh * J) * focal;

  if (focal == -1 && neigh == -4) {
    en_diff = en_diff + inertia;
  }

  return en_diff;
}

// [[Rcpp::export]]
NumericMatrix flip_glauber_rcpp(const NumericMatrix& input_matrix, const double& B, const double& J, const IntegerVector& rxs, const IntegerVector& rys, const NumericVector& rns, const int& n_rows, const int& n_cols, const double& inertia) {

  NumericMatrix output_matrix = clone(input_matrix);
  int n = rns.size();

  for (int i = 0; i < n; i++) {
    int rx = rxs[i], ry = rys[i];
    double rn = rns[i];

    double nb = output_matrix(mod(rx, n_rows), ry - 1) +
      output_matrix(mod((rx - 2), n_rows), ry - 1) +
      output_matrix(rx - 1, mod(ry, n_cols)) +
      output_matrix(rx - 1, mod((ry - 2), n_cols));

    double fo = output_matrix(rx - 1, ry - 1);
    double en_diff = energy_diff3(fo, nb, B, J, inertia);
    double P = 1 / (1 + exp(en_diff));

    if (P > rn) {
      output_matrix(rx - 1, ry - 1) = -fo;
    }
  }

  return output_matrix;
}

// [[Rcpp::export]]
NumericMatrix flip_metropolis2_rcpp(const NumericMatrix& input_matrix, const double& B, const double& J, const IntegerVector& rxs, const IntegerVector& rys, const NumericVector& rns, const int& n_rows, const int& n_cols, const double& inertia) {

  NumericMatrix output_matrix = clone(input_matrix);
  int n = rns.size();

  for (int i = 0; i < n; i++) {
    int rx = rxs[i], ry = rys[i];
    double rn = rns[i];

    double nb = output_matrix(mod(rx, n_rows), ry - 1) +
      output_matrix(mod((rx - 2), n_rows), ry - 1) +
      output_matrix(rx - 1, mod(ry, n_cols)) +
      output_matrix(rx - 1, mod((ry - 2), n_cols));

    double fo = output_matrix(rx - 1, ry - 1);
    double en_diff = energy_diff3(fo, nb, B, J, inertia);

    if (en_diff <= 0){
      output_matrix(rx - 1, ry - 1) = -fo;
    } else {
      double p = exp(-en_diff);
      if (rn < p){
        output_matrix(rx - 1, ry - 1) = -fo;
      }
    }
  }

  return output_matrix;
}
