#include "isotopesanCliqueMSR.h"

// [[Rcpp::export]]
Rcpp::DataFrame returnIsotopes(Rcpp::DataFrame dfclique, int maxCharge, double ppm, double isom = 1.003355) {
  isoData isoClique = readisoData(dfclique);
  isoReturn result = getIsotopes(isoClique, maxCharge, ppm, isom);
  return Rcpp::DataFrame::create(Rcpp::Named("pfeature") = result.pfeature,
				 Rcpp::Named("ifeature") = result.ifeature,
				 Rcpp::Named("pcharge") = result.pcharge,
				 Rcpp::Named("icharge") = result.icharge);
}
