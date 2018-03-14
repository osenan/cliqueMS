#ifndef ISOTOPESANCLIQUEMS_H_
#define ISOTOPESANCLIQUEMS_H_

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>


class isoData {
 public:
  std::vector<double> mz;
  std::vector<int> feature;
};
  
bool errorRange (double mz1, double mz2, double reference, double ppm)
{
  bool result = false;
  long double error = (std::abs(mz2 -mz1 -reference))/(mz1+reference);
  if( error <= sqrt(2.0)*ppm*0.000001 )
    result = true;
  return result;
}

isoData readisoData(Rcpp::DataFrame dfclique)
{
  isoData isoD;
  Rcpp::NumericVector vmz = dfclique["mz"];
  Rcpp::NumericVector vfeature = dfclique["feature"];
  for(int index = 0; index < vmz.size(); index++) {
    isoD.feature.push_back(vfeature[index]);
    isoD.mz.push_back(vmz[index]);
  }
  return isoD;
}

class isoTest {
public:
  bool isIso = false;
  int pcharge;
  int icharge;
};

isoTest isIsotope(double mz1, double mz2, int maxCharge, double ppm, double isom = 1.003355)
{
  double cmz1, cmz2;
  int charge1;
  int charge2;
  isoTest currentIso, finalIso;
  for(charge1 = 1; charge1 <= maxCharge; charge1++) {
    for(charge2 = 1; charge2 <= maxCharge; charge2++) {
      cmz1 = mz1*charge1;
      cmz2 = mz2*charge2;
      // if isotope mass is bigger than parental mass
      if(cmz2 > cmz1) {
	currentIso.isIso = errorRange(cmz1, cmz2, isom, ppm);
	if( currentIso.isIso == true) {
	  //std::cout << cmz1 << " " << cmz2 << " " << "\n";
	  finalIso.isIso = true;
	  finalIso.pcharge = charge1;
	  finalIso.icharge = charge2;
	}
      }
    }
  }
  return finalIso;
}

class isoReturn {
 public:
  std::vector <int> pfeature;
  std::vector <int> ifeature;
  std::vector <int> pcharge;
  std::vector <int> icharge;
};

isoReturn getIsotopes(isoData isoD, int maxCharge, double ppm, double isom = 1.003355)
{
  isoReturn result;
  double mz1, mz2;
  isoTest iTest;
  for(unsigned int id1 = 0; id1 < isoD.mz.size() ; id1++) {
    for(unsigned int id2 = 1; id2 < isoD.mz.size() ; id2++) {
      if(id2 > id1) {
	mz1 = isoD.mz[id1];
	mz2 = isoD.mz[id2];
	iTest = isIsotope(mz1, mz2, maxCharge, ppm, isom);
	if( iTest.isIso == true ) {
	  //std::cout << id1 << " " << id2 << "\n";
	  result.pfeature.push_back(isoD.feature[id1]);
	  result.ifeature.push_back(isoD.feature[id2]);
	  result.pcharge.push_back(iTest.pcharge);
	  result.icharge.push_back(iTest.icharge);
	}
      }
    }
  }
  return result;
}
#endif

