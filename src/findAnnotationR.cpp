#include "annotationCliqueMSR.h"

// [[Rcpp::export]]
Rcpp::DataFrame returnAnnotation(Rcpp::DataFrame dfclique, Rcpp::DataFrame dfaddlist, int topmassf = 1, int topmasstotal = 10, unsigned int sizeanG = 20, double tol = 0.00001, double filter = 0.0001, double emptyS = -6) {

  double newmass = -10.0;
  outputAn outAn = getAnnotation(dfclique, dfaddlist, topmassf, topmasstotal, sizeanG, tol, filter, emptyS, newmass);
  std::vector<double> vmass1, vmass2, vmass3, vmass4, vmass5, vscore1, vscore2, vscore3, vscore4, vscore5;
  std::vector<std::string> van1, van2, van3, van4, van5;

  for(std::vector<int>::iterator itv = outAn.features.begin(); itv != outAn.features.end(); itv++) {
    vmass1.push_back(outAn.mass1[*itv]);
    vmass2.push_back(outAn.mass2[*itv]);
    vmass3.push_back(outAn.mass3[*itv]);
    vmass4.push_back(outAn.mass4[*itv]);
    vmass5.push_back(outAn.mass5[*itv]);
    vscore1.push_back(outAn.score1[*itv]);
    vscore2.push_back(outAn.score2[*itv]);
    vscore3.push_back(outAn.score3[*itv]);
    vscore4.push_back(outAn.score4[*itv]);
    vscore5.push_back(outAn.score5[*itv]);
    van1.push_back(outAn.an1[*itv]);
    van2.push_back(outAn.an2[*itv]);
    van3.push_back(outAn.an3[*itv]);
    van4.push_back(outAn.an4[*itv]);
    van5.push_back(outAn.an5[*itv]);
  }
  return Rcpp::DataFrame::create(
				 Rcpp::Named("feature") = outAn.features,

				 Rcpp::Named("mass1") = vmass1,
				 Rcpp::Named("an1") = van1,
				 Rcpp::Named("score1") = vscore1,

				 Rcpp::Named("mass2") = vmass2,
				 Rcpp::Named("an2") = van2,
				 Rcpp::Named("score2") = vscore2,

				 Rcpp::Named("mass3") = vmass3,
				 Rcpp::Named("an3") = van3,
				 Rcpp::Named("score3") = vscore3,

				 Rcpp::Named("mass4") = vmass4,
				 Rcpp::Named("an4") = van4,
				 Rcpp::Named("score4") = vscore4,

				 Rcpp::Named("mass5") = vmass5,
				 Rcpp::Named("an5") = van5,
				 Rcpp::Named("score5") = vscore5
				 );
}
