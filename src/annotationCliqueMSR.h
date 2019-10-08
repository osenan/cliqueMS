// [[Rcpp::depends(BH)]]
#ifndef ANNOTATIONCLIQUEMS_H_
#define ANNOTATIONCLIQUEMS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <iomanip>
#include <Rcpp.h>

#include <boost/functional/hash.hpp>



typedef boost::hash<std::pair<double, double> > Hashd;
typedef std::unordered_set<std::pair<double, double>, Hashd> pairM;
typedef boost::hash<std::pair<double, double> > Hashi;
typedef std::unordered_set<std::pair<int, int>, Hashi> pairA;

class annotDF {
 public:
  std::vector<double> mz;
  std::vector<int> features;
  std::vector<int> charge;
};

annotDF readDF(Rcpp:: DataFrame dfclique)
{
  annotDF annotdf;
  // read data.frame from R
  Rcpp::NumericVector vmz = dfclique["mz"];
  Rcpp::NumericVector vfeature = dfclique["feature"];
  Rcpp::NumericVector vcharge = dfclique["charge"];
  // copy data.frame data to c++ annotDF class
  for(int index = 0; index < vmz.size(); index++) {
    annotdf.mz.push_back(vmz[index]);
    annotdf.features.push_back(vfeature[index]);
    annotdf.charge.push_back(vcharge[index]);
  }
  return annotdf;
}

class adInfo {
 public:
  double freq;
  double massDiff;
  int numMol;
  int charge;
};


class rawadList {
 public:
  std::unordered_map <std::string, adInfo> rawlist;
  std::vector <std::string> addorder;
};
  
rawadList readrawList(Rcpp:: DataFrame dfaddlist)
{
  rawadList rawL;
  adInfo adI;
  std::string add;
  Rcpp::StringVector vadd = dfaddlist["adduct"];
  Rcpp::NumericVector vlog10freq = dfaddlist["log10freq"];
  Rcpp::NumericVector vmassdiff = dfaddlist["massdiff"];
  Rcpp::NumericVector vnummol = dfaddlist["nummol"];
  Rcpp::NumericVector vcharge = dfaddlist["charge"];
  for(int index = 0; index < vadd.size(); index++) {
    adI.freq = vlog10freq[index];
    adI.massDiff = vmassdiff[index];
    adI.numMol = vnummol[index];
    adI.charge = vcharge[index];
    add = vadd[index];
    rawL.rawlist[add] = adI; // store adduct with all its information about freq, charge, massdiff and number of molecules
    rawL.addorder.push_back(add);
  }
  return rawL;
}

std::vector<double> getScoreAddlist (rawadList& rList) {
  // function to sort the list of adducts by score
  std::vector <double> vScore;
  for(std::unordered_map<std::string, adInfo>::iterator itl = rList.rawlist.begin(); itl != rList.rawlist.end(); itl++)
    vScore.push_back(itl->second.freq);
  std::sort(vScore.begin(), vScore.end());
  return(vScore);
}

class annotData {
 public:
  std::unordered_map <int, int> features;
  std::unordered_map <double, std::vector<std::pair<int,std::string>> > massList;
  std::unordered_map <int, std::vector<double> > feat2mass;
  std::unordered_map <int, std::vector<int> > anGroups;
  std::unordered_map <int, std::unordered_map <double, std::vector<std::pair<int,std::string>> > > anGroup2mass;
};

std::unordered_map <int, std::string> getAlladducts(double mass, double tol, int idn, adInfo currentAdd, annotDF& mzdf, rawadList rList) {
  unsigned int idnmass = idn; //index to start the search in the row of adducts
  int isoCharge; // charge set to isotopic features
  double mapmassDiff, mzDiff, error, lowerbound, upperbound;
  adInfo adI, adItest;
  std::string addTest;
  std::map <double,std::string> massMap;
  std::unordered_map <int, std::string> adduMap;
  for(std::vector <std::string>::iterator ita = rList.addorder.begin(); ita!= rList.addorder.end() ; ita++) {
    adI = rList.rawlist[*ita];
    mapmassDiff = -mass + (mass*adI.numMol + adI.massDiff)/std::abs(adI.charge);
    massMap[mapmassDiff] = *ita;
  }
  std::map<double,std::string>::iterator lowerboundp = massMap.begin();
  lowerbound = lowerboundp->first -lowerboundp->first*0.10;
  std::map<double,std::string>::reverse_iterator upperboundp = massMap.rbegin();
  upperbound = upperboundp->first +upperboundp->first*0.10;
  // first see if there is any previous row prior to the one that we start the search
  while( (mzdf.mz[idnmass]-mass) < lowerbound) {
    if( idnmass == 0 )
      break;
    idnmass --;
  }
  // search for all adducts of the mass in the df
  for(unsigned int idnloop = idnmass; idnloop < mzdf.mz.size() ; idnloop++ ) {
    mzDiff = mzdf.mz[idnloop] -mass;
    isoCharge = mzdf.charge[idnloop];
  // if massdifference is bigger than the largest mass difference in the adduct list
  // it is not possible to find more adducts
    std::map<double,std::string>::iterator itm;
      for( itm = lowerboundp; itm != massMap.end() ; itm++ ) {
	if(isoCharge != 0) {
	  addTest = itm->second;
	  adItest = rList.rawlist[addTest];
	  if(isoCharge == std::abs(adItest.charge)) {
	    // if testing add charge is equal to feature charge
	    error = std::abs(mzDiff - itm->first)/mass;
	  } else {
	    // if not, this annotation is not possible
	    error = 10*tol*sqrt(2);
	  }
	} else {
	  // if there is no charge set, proceed as normal
	  error = std::abs(mzDiff - itm->first)/mass;
	}
    // if the error is smaller than the tolerance, accept that adduct
	if( error < tol*sqrt(2) )
	  adduMap[mzdf.features[idnloop]] = itm->second;
      }
      // move lowerbound if the mass difference is larger than the lowerbound
      if( mzDiff > lowerboundp->first )
	std::advance(lowerboundp,1);
      if( mzDiff > upperbound )
	break;
  }
  return adduMap;
}

class Component {
 public:
  std::unordered_set<double> mass;
  std::unordered_set<int> feature;
};

void getComponent(std::unordered_set<double>& setm, std::unordered_set<int>& extraf, annotData& annotD, Component& comp) {
  std::unordered_set<double> extram;
  std::unordered_set<int> newf;
  for(std::unordered_set<int>::iterator itf = extraf.begin(); itf!= extraf.end(); itf++) {
    comp.feature.insert(*itf); // insert this feature if it is not in the component;
    std::unordered_set<int>::iterator it;
    for(std::vector<double>::iterator itm = annotD.feat2mass[*itf].begin(); itm!= annotD.feat2mass[*itf].end(); itm++) {
      // if in this feature there is mass on the setm
      if(setm.find(*itm) != setm.end()) {
	// include the mass in the component if it is not yet
	auto mnew = comp.mass.insert(*itm);
	// if this mass is new, search for new features
	if(mnew.second == true) {
	  for(std::vector<std::pair<int,std::string>>::iterator itvm = annotD.massList[*itm].begin();
	      itvm != annotD.massList[*itm].end(); itvm++) {
	    // if some feature is not in the component feature set, include it in the new features
	    if(comp.feature.find(itvm->first) == comp.feature.end())
	      newf.insert(itvm->first);
	  }
	}
      }
    }
  }
  // if there are no new features to check end
  if(newf.size() > 0)
    getComponent(setm, newf, annotD, comp);
}

void getComponentanG(std::unordered_set<int>& extraf, annotData& annotD, Component& comp) {
  std::unordered_set<int> newf;
  for(std::unordered_set<int>::iterator itf = extraf.begin(); itf!= extraf.end(); itf++) {
    comp.feature.insert(*itf); // insert this feature if it is not in the component;
    for(std::vector<double>::iterator itm = annotD.feat2mass[*itf].begin(); itm!= annotD.feat2mass[*itf].end(); itm++) {
      // if in this feature there is mass on the setm
      // include the mass in the component if it is not yet
      auto mnew = comp.mass.insert(*itm);
      // if this mass is new, search for new features
      if(mnew.second == true) {
	for(std::vector<std::pair<int,std::string>>::iterator itvm = annotD.massList[*itm].begin();
	    itvm != annotD.massList[*itm].end(); itvm++) {
	  // if some feature is not in the component feature set, include it in the new features
	  if(comp.feature.find(itvm->first) == comp.feature.end())
	    newf.insert(itvm->first);
	}
      }
    }
  }

  // if there are no new features to check end
  if(newf.size() > 0)
    getComponentanG(newf, annotD, comp);
}


std::unordered_map<int, Component> getanGcomp(annotData& annotD) {
  std::unordered_map<int, Component> mapC;
  int id = 0;
  std::unordered_set<int> setf;
  // 1 - include all the features in this clique
  for(std::unordered_map<int, int>::iterator itv = annotD.features.begin(); itv != annotD.features.end(); itv++)
    setf.insert(itv->first);
  // 2 - Compute wich anGroups exist
  while(setf.size() != 0) {
    Component comp;
    std::unordered_set<int> extraf;
    std::unordered_set<int>::iterator it = setf.begin();
    extraf.insert(*it);
    getComponentanG(extraf, annotD, comp);
    // drop features if comp created
    if(comp.feature.size() > 0) {
      for(std::unordered_set<int>::iterator itf = comp.feature.begin(); itf != comp.feature.end(); itf++)
	setf.erase(*itf);
      mapC[id] = comp;
      id++;
    }
  }
  return mapC;
}

bool compareMasses(annotData& annotD, double m1, double m2, int anGroup) {
  //m1 always is the mass to drop
  bool result = false;
  unsigned int count = 0;
  for(std::vector<std::pair<int, std::string> >::iterator it1 = annotD.anGroup2mass[anGroup][m1].begin();
      it1 != annotD.anGroup2mass[anGroup][m1].end(); it1++) {
    for(std::vector<std::pair<int, std::string> >::iterator it2 = annotD.anGroup2mass[anGroup][m2].begin();
	it2 != annotD.anGroup2mass[anGroup][m2].end(); it2++) {
      if( *it1 == *it2 ) 
	count++;
    }
  }
  if( count == annotD.anGroup2mass[anGroup][m1].size() )
    result = true;
  return result;
}

std::unordered_set<double> getRepMasses(int anGroup, annotData& annotD, double filter = 0.0001) {
  double error, m1, m2;
  // pair of masses to check, they are ordered from small to big value
  std::pair<double,double> pairmass; // the second is the bigger
  pairM badmPair;
  std::unordered_set<double> badmasses;
  for(std::unordered_map<double, std::vector<std::pair<int,std::string>> >::iterator itm1 = annotD.anGroup2mass[anGroup].begin();
      itm1 != annotD.anGroup2mass[anGroup].end(); itm1++) {
    for(std::unordered_map<double, std::vector<std::pair<int,std::string>> >::iterator itm2 = annotD.anGroup2mass[anGroup].begin();
      itm2 != annotD.anGroup2mass[anGroup].end(); itm2++) {
      if(itm1 != itm2) {
	// now check if masses are very similar
	m1 = itm1->first;
	m2 = itm2->first;
	error = std::abs(m1-m2)/m1;
	// if the error is smaller than the filter, we have to check the adducts
	if( error < filter*sqrt(2) ) {
	  if(m2 > m1) {
	    pairmass.first = m1;
	    pairmass.second = m2;
	  } else {
	    pairmass.first = m1;
	    pairmass.second = m2;
	  }
	  badmPair.insert(pairmass);
	}
      }
    }
  }
  // now filter the masses that are similar
  for(pairM::iterator itmp = badmPair.begin(); itmp != badmPair.end(); itmp ++) {
    int size1, size2;
    bool comp = false;
    size1 = annotD.anGroup2mass[anGroup][itmp->first].size();
    size2 = annotD.anGroup2mass[anGroup][itmp->second].size();
    m1 = itmp->first;
    m2 = itmp->second;
    if(size1 == size2) {
      // if masses are equal, we eventually drop m1
      comp = compareMasses(annotD, m1, m2, anGroup);
      if(comp == true) 
	badmasses.insert(m1);
    } else {
      if(size1 > size2) {
	// size1 is bigger, drop m2 if conditions met
	comp = compareMasses(annotD, m2, m1, anGroup);
	if(comp == true) 
	  badmasses.insert(m2);
      } else {
	//size 2 is bigger, drop m1 if conditions met
	comp = compareMasses(annotD, m1, m2, anGroup);
	if(comp == true) 
	  badmasses.insert(m1);
      }
    }
  }
  return(badmasses);
}


void createanGroup2mass(annotData& annotD, std::unordered_map<int, Component>& anGcomp, double filter = 0.0001) {
  // first create anGroup2mass object
  for(std::unordered_map<int, Component>::iterator itc = anGcomp.begin(); itc != anGcomp.end(); itc++) {
    for(std::unordered_set<double>::iterator its = anGcomp[itc->first].mass.begin();
	its != anGcomp[itc->first].mass.end(); its++)
      annotD.anGroup2mass[itc->first][*its] = annotD.massList[*its];
  }

  
  // second filter masses with very similar mass and the same adducts, or less adducts in one case
  std::unordered_set<double> setm;
  for(std::unordered_map <int, std::unordered_map<double, std::vector<std::pair<int,std::string>> > >::iterator itg =
	annotD.anGroup2mass.begin(); itg!= annotD.anGroup2mass.end(); itg++) {
    setm = getRepMasses(itg->first, annotD, filter);
    // if there are repeated masses drop these masses
    if(setm.size() > 0) {
      for(std::unordered_set<double>::iterator itms = setm.begin(); itms != setm.end(); itms++) {
	annotD.anGroup2mass[itg->first].erase(*itms);
	annotD.massList.erase(*itms); // also drop them in the massList
	for(std::vector<int>::iterator itf = annotD.anGroups[itg->first].begin();
	    itf != annotD.anGroups[itg->first].end(); itf++) {
	  auto mf = std::find(annotD.feat2mass[*itf].begin(), annotD.feat2mass[*itf].end(), *itms);
	  if(mf != annotD.feat2mass[*itf].end() )
	    annotD.feat2mass[*itf].erase(mf); // and drop them from the mass assigned to each feature
	}
      }
    }
  }
}
      
annotData getannotData(rawadList rList, annotDF& mzdf, double tol = 0.00001, double filter = 0.0001)
{
  adInfo currentAdd;
  annotData annotD;
  unsigned int idn;
  // initialize annotD
  for(idn = 0; idn < mzdf.features.size(); idn++)
    annotD.features[mzdf.features[idn]] = -1;
  
  std::pair<int, std::string> massPair;
  double mz, mass;
  int isocharge;
  std::unordered_map<int,std::string> adduMap;
  for(idn = 0; idn < (mzdf.mz.size() -1); idn++) {
    mz = mzdf.mz[idn]; // assign m/z value for current idn position
    isocharge = mzdf.charge[idn]; // assign charge value, only available for isotopic features
    for(std::vector<std::string>::iterator ita = rList.addorder.begin(); ita != rList.addorder.end(); ita++) {
      currentAdd = rList.rawlist[*ita];
      if(isocharge != 0) {
        // if mz is an isotopic feature with assigned charge
	mass = -1.0;
	if(std::abs(currentAdd.charge) == isocharge) {
	  // only annotate adduct if the charge of the isotopic feature is equal to the charge of the putative adduct
	  mass = (mz*std::abs(currentAdd.charge) - currentAdd.massDiff)/currentAdd.numMol;
	}
      } else {
	mass = (mz*std::abs(currentAdd.charge) - currentAdd.massDiff)/currentAdd.numMol;
      }
      if(mass > 0) {
	mass = round(mass*10000)/10000;
	if(annotD.massList.find(mass) == annotD.massList.end() ) {
	  // if this mass is not on the mass list, search for all the adducts of that mass
	  adduMap = getAlladducts(mass, tol, idn, currentAdd, mzdf, rList);
	  // if there is more than one adduct:
	  if( adduMap.size() > 1 ) {
	    for(std::unordered_map<int, std::string>::iterator itu = adduMap.begin();
		itu != adduMap.end(); itu++) {
	      // update massList
	      massPair.first = itu->first;
	      massPair.second = itu->second;
	      annotD.massList[mass].push_back(massPair);
	      // update at which features is this mass
	      annotD.feat2mass[itu->first].push_back(mass);
	    }
	  }
	}
      }
    }
  }
  // create anGroups
  std::unordered_map<int, Component> anGcomp =  getanGcomp(annotD);
  for(std::unordered_map<int, Component>::iterator itc = anGcomp.begin(); itc != anGcomp.end(); itc++) {
    for(std::unordered_set<int>::iterator itf = anGcomp[itc->first].feature.begin();
	itf != anGcomp[itc->first].feature.end(); itf++) {
      annotD.anGroups[itc->first].push_back(*itf);
      annotD.features[*itf] = itc->first;
    }
  }
  // create and filter angroup2mass
  createanGroup2mass(annotD, anGcomp, filter);
  return annotD;
}

bool compare(const std::pair<double, double>& pair1, const std::pair<double, double>& pair2) {
  // sort a vector on pairs in descending orders, in that case sort score and mass pair
  return pair1.first > pair2.first;
}

bool compareint(const std::pair<double, int>& pair1, const std::pair<double, int>& pair2) {
  // sort a vector on pairs in descending orders, in that case sort score and mass pair
  return pair1.first > pair2.first;
}


std::vector< std::pair<double,double> > sortMass (annotData& annotD, int feature,
						  std::unordered_map< double, std::pair<double, double> > mass2score, int n) {
  std::vector< std::pair<double,double> > allM;
  std::vector< std::pair<double, double> > topV;
  for(std::vector<double>::iterator itm = annotD.feat2mass[feature].begin();
      itm != annotD.feat2mass[feature].end(); itm++)
    allM.push_back(mass2score[*itm]);
  // sort mass vector according to score
  sort(allM.begin(), allM.end(), compare);
  // select the top "n" masses
  for(int id = 0; id < n; id++) {
    if(id < allM.size()) // not add more masses in case that for that feature are less than "n" top masses
      topV.push_back(allM[id]);
  }
  return topV;
}

std::unordered_set<double> getTopScoringMasses(annotData& annotD, int anG, rawadList rList, int nfeature = 1,
					       int ntotal = 10, double emptyS = -6) {
  // we want to get the "n" top scoring masses for each feature of the annotationGroup 
  std::unordered_set<double> setm;
  std::vector<std::pair<double, double> > topF, allM, topT;
  std::unordered_map<double, std::pair<double, double> > mass2score;
  std::pair<double, double> mass2scoreEntry;
  adInfo adI;
  //double minsert;
  // First compute the score for all the masses in the anGroup
  for(std::unordered_map<double, std::vector<std::pair<int,std::string>> >::iterator itm = annotD.anGroup2mass[anG].begin();
      itm != annotD.anGroup2mass[anG].end(); itm++) {
      double score = 0;
      for(std::vector<std::pair<int, std::string>>::iterator ita = annotD.anGroup2mass[anG][itm->first].begin();
	  ita != annotD.anGroup2mass[anG][itm->first].end(); ita++) {
	adI = rList.rawlist[ita->second];
	score += adI.freq; // compute the score according to the log frequence of each adduct in the adduct list
      }
      //add the compensation for empty annotations
      score += emptyS*(annotD.anGroups[anG].size() - annotD.anGroup2mass[anG][itm->first].size());
      mass2scoreEntry.first = score;
      mass2scoreEntry.second = itm->first;
      mass2score[itm->first] = mass2scoreEntry;
  }
  // Second select the top masses independently of the feature for that group
  for(std::unordered_map<double, std::pair<double, double> >::iterator itm1 = mass2score.begin();
      itm1 != mass2score.end(); itm1++)
    allM.push_back(itm1->second);
  sort(allM.begin(), allM.end(), compare); // sort this masses
  for(int id = 0; id < ntotal; id++) {
    if(id < allM.size()) // not add more masses in case that there are less than "n" top masses
      topT.push_back(allM[id]);
  }
  // Add this masses to the total set of masses
  for(std::vector<std::pair<double,double>>:: iterator itv = topT.begin(); itv != topT.end(); itv++)
    setm.insert(itv->second);
  // Third, for each feature, select the "n" top scoring masses
  for(std::vector<int>::iterator itf = annotD.anGroups[anG].begin(); itf != annotD.anGroups[anG].end(); itf++) {
    topF = sortMass(annotD, *itf, mass2score, nfeature);
    for(std::vector<std::pair<double,double>>:: iterator itv = topF.begin(); itv != topF.end(); itv++)
      setm.insert(itv->second);
  }
  return setm;
}


class Annotation {
 public:
  double score = 0;
  std::unordered_map<int, std::pair<double, std::string> > annotation;
};

Annotation annotateMass(annotData& annotD, std::unordered_set<int> features, rawadList rList,
						     std::unordered_set<double> setm, double emptyS = -6) {
  Annotation an, anfree;
  anfree.score = 0;
  an.score = 0;
  int count;
  double score, topmass;
  adInfo adI;
  std::vector<std::pair<double, double> > mass2score; // first is the score, second is the mass
  std::pair<double, double> mass2scoreEntry; // entry for mass2score vector
  std::pair<double, std::string> anEntry; // entry for an
  std::unordered_set<int> freef;
  // 1 - for each mass in setm compute the score for the features if that features are more than in two positions
  for(std::unordered_set<double>::iterator itm = setm.begin(); itm != setm.end(); itm++) {
    score = 0;
    count = 0;
    for(std::unordered_set<int>::iterator itf = features.begin(); itf != features.end(); itf++) {
      //check if this feature is in the massList
      for(std::vector<std::pair<int, std::string> >::iterator itp = annotD.massList[*itm].begin();
	  itp != annotD.massList[*itm].end(); itp++) {
	// if this mass contains the feature itf, compute the score
	if(itp->first == *itf) {
	  adI = rList.rawlist[itp->second];
	  score += adI.freq;
	  count++;
	}
      }
    }
    if(count > 1) {
      mass2scoreEntry.first = score;
      mass2scoreEntry.second = *itm;
      mass2score.push_back(mass2scoreEntry);
    }
  }
  
  // if there is no annotation for the features, return an empty score
  if( mass2score.size() < 1 ) {
    an.score = emptyS*features.size();
    for(std::unordered_set<int>::iterator itf = features.begin(); itf != features.end(); itf++) {
      anEntry.first = 0;
      anEntry.second = "";
      an.annotation[*itf] = anEntry;
    }
    return(an);
  }
  // 2 - sort annotation and select the adducts of that annotation
  sort(mass2score.begin(), mass2score.end(), compare);
  topmass = mass2score[0].second;
  an.score = mass2score[0].first;
  an.score += -10; // add the log compensation of -10 for each new mass;

  
  // search again for adduct in that feature and include in annotation
  for(std::unordered_set<int>::iterator itf = features.begin(); itf != features.end(); itf++) {
    for(std::vector<std::pair<int, std::string> >::iterator itp = annotD.massList[topmass].begin();
	itp != annotD.massList[topmass].end(); itp++) {
      // if this mass contains the feature itf, include it in the annotation
      if(itp->first == *itf) {
	anEntry.first = topmass;
	anEntry.second = itp->second;
      	an.annotation[*itf] = anEntry; // mass and adduct
      }
    }
  }
  // now put the features out of the annotation in the freef
  std::unordered_map<int, std::pair<double, std::string> >::iterator ita;
  for(std::unordered_set<int>::iterator itf = features.begin(); itf != features.end(); itf++) {
    ita = an.annotation.find(*itf);
    if( ita == an.annotation.end() )
      freef.insert(*itf);
  }
      
  if( freef.size() > 1 ) {
    setm.erase(topmass);
    anfree = annotateMass(annotD, freef, rList, setm, emptyS);
    for(std::unordered_map<int, std::pair<double, std::string> >::iterator itan = anfree.annotation.begin();
	itan != anfree.annotation.end(); itan++) {
      // add the additional score and annotation after the recursive call
      an.annotation[itan->first] = itan->second;
    }
    an.score += anfree.score;
  } else {
    if( freef.size() == 1 ) {
      an.score += emptyS;
      auto it = freef.begin();
      anEntry.first = 0;
      anEntry.second = "";
      an.annotation[*it] = anEntry;
    }
  }
  return(an);
}

bool compareAnnotation(int id1, int id2, std::unordered_map<int, Annotation>& annotations) {
  bool result = false;
  int count = annotations[id1].annotation.size();
  for(std::unordered_map<int, std::pair<double, std::string>>::iterator itf = annotations[id1].annotation.begin();
      itf != annotations[id1].annotation.end(); itf++) {
    // one less for the count if the annotation for that feature is the same
    if(itf->second == annotations[id2].annotation[itf->first])
      count--;
  }
  // if all annotations are the same, drop than annotation
  if(count == 0)
    result = true;
  return result;
}
      

void dropRepeatedAnnotations(std::unordered_map<int, Annotation>& annotations) {
  pairA badA;
  std::pair<int,int> badAEntry;
  double score1, score2;
  int id1,id2;
  std::unordered_set<int> dropA;
  for(std::unordered_map<int, Annotation>::iterator itm1 = annotations.begin(); itm1 != annotations.end(); itm1++) {
    score1 = itm1->second.score;
    for(std::unordered_map<int, Annotation>::iterator itm2 = annotations.begin(); itm2 != annotations.end(); itm2++) {
      // if is not the same annotation
      if(itm1->first != itm2->first) {
	score2 = itm2->second.score;
	if(score1 == score2) {
	  // if itm1->first annotation is bigger
	  if(itm1->first > itm2->first) {
	    badAEntry.first = itm2->first;
	    badAEntry.second = itm1->first;
	    badA.insert(badAEntry);
	  } else{
	    badAEntry.first = itm1->first;
	    badAEntry.second = itm2->first;
	    badA.insert(badAEntry);
	  }
	}
      }
    }
  }
  for( pairA::iterator itp = badA.begin(); itp != badA.end(); itp++ ) {
    id1 = itp->first;
    id2 = itp->second;
    bool compareA = compareAnnotation(id1, id2, annotations);
    // if the two annotations are exactly the same, drop the annotation with a bigger id
    if(compareA == true)
      dropA.insert(id2);
  }
  // if there are equal annotations, drop them from annotations map
  if(dropA.size() > 0) {
    for(std::unordered_set<int>::iterator its = dropA.begin(); its != dropA.end(); its++)
      annotations.erase(*its);
  }
}

std::unordered_map<int, Annotation> AnnotateMassLoop(annotData& annotD, std::unordered_set<int> features, rawadList rList,
						     std::unordered_set<double> setm, double emptyS = -6) {
  std::unordered_map<int, Annotation> annotations;
  int id = 0;
  adInfo adI;
  // 1 - For each mass get the annotations that coincide with features
  for(std::unordered_set<double>::iterator itm = setm.begin(); itm != setm.end(); itm++) {
    std::unordered_set<int> freef;
    std::unordered_set<double> setm2 = setm;
    Annotation an, anfree;
    an.score = 0;
    anfree.score = 0;
    std::pair<double, std::string> anEntry; // entry for an
    for(std::unordered_set<int>::iterator itf = features.begin(); itf != features.end(); itf++) {
      for(std::vector<std::pair<int, std::string> >::iterator itp = annotD.massList[*itm].begin();
	itp != annotD.massList[*itm].end(); itp++) {
	// if this mass contains the feature itf, include it in the annotation
    	if(itp->first == *itf) {
	  anEntry.first = *itm;
	  anEntry.second = itp->second;
	  an.annotation[*itf] = anEntry; // mass and adduct
	  adI = rList.rawlist[itp->second];
	  an.score += adI.freq;
	}
      }
    }

    // 2 - Put the features out of the annotation of this mass *itm in the freef
    std::unordered_map<int, std::pair<double, std::string> >::iterator ita;
    for(std::unordered_set<int>::iterator itf = features.begin(); itf != features.end(); itf++) {
      ita = an.annotation.find(*itf);
      if( ita == an.annotation.end() )
	freef.insert(*itf);
    }
    
  // 3 - Execute AnnotateMass with the freef features
    if( freef.size() > 1 ) {
      setm2.erase(*itm);
      anfree = annotateMass(annotD, freef, rList, setm2, emptyS);
      for(std::unordered_map<int, std::pair<double, std::string> >::iterator itan = anfree.annotation.begin();
	  itan != anfree.annotation.end(); itan++) {
      // add the additional score and annotation after the recursive call
	an.annotation[itan->first] = itan->second;
      }
      an.score += anfree.score;
    } else {
      if( freef.size() == 1 ) {
	an.score += emptyS;
	auto it = freef.begin();
	anEntry.first = 0;
	anEntry.second = "";
	an.annotation[*it] = anEntry;
      }
    }
  // 4 - Include this annotation in the map
    annotations[id] = an;
    id++;
  }
  // 5 - Filter equal annotations

  dropRepeatedAnnotations(annotations);
  return annotations;
}

std::vector<int> sortAnnotations(std::unordered_map<int, Annotation>& annotations, int top) {
  std::vector<std::pair<double, int> > allAn;
  std::pair<double, int> allAnEntry;
  std::vector<int> topAn;
  
  for(std::unordered_map<int, Annotation>::iterator ita = annotations.begin(); ita != annotations.end(); ita++) {
    allAnEntry.first = annotations[ita->first].score;
    allAnEntry.second = ita->first;
    allAn.push_back(allAnEntry);
  }
  sort(allAn.begin(), allAn.end(), compareint);
  for(int id = 0; id < top; id++) {
    if(id < allAn.size() )
      topAn.push_back(allAn[id].second);
  }
  return topAn;
}

double computeMaxScore(std::vector<double>& vScore, int annotsize, double newmass = -10.0) {
  double score = 0;
  double completeroundscore,remainderroundscore = 0;
  int completeround = annotsize/vScore.size();
  int remainderround = annotsize%vScore.size();
  std::vector<double>::reverse_iterator ritv;
  // compute score by number of complete rounds
  for(ritv = vScore.rbegin(); ritv != vScore.rend(); ritv++)
    completeroundscore += *ritv;
  // compute score by adducts on the remainder
  ritv = vScore.rbegin();
  for(int i = 0; i < remainderround; i++) {
    remainderroundscore += *ritv;
    ritv++;
  }
  //the final score is the number of loops with the total list, plus the number of extramasses, and the remainder adduct not complete
  score = (completeround*completeroundscore) + remainderroundscore + (completeround*newmass);
  return(score);
}

void normalizeAnnotation(Annotation& an, std::vector<double>& vScore, double newmass = -10.0,
			 double emptyS = -6, int newmassSize = 8.0) {
  double maxscore, minscore, oldscore, newscore = 0;
  oldscore = an.score;
  maxscore = computeMaxScore(vScore, an.annotation.size(), newmass);
  // computation of min score is as if all the annotation is empty score
  // plus a compensation of a new mass each 8 features
  minscore = an.annotation.size()*emptyS + newmass*(an.annotation.size()/newmassSize);
  newscore = 100*(oldscore - minscore)/(maxscore - minscore); // taken from the linear interpolation formula
  if(newscore < 0) // there are cases where new score is below 0, in those cases the normalized score should be scaled to zero
    newscore = 0;
  an.score = round(10000*(newscore))/10000;
}

std::unordered_map<int, Component> getSeparateComp(std::unordered_set<double> setm, annotData& annotD, int anG) {
  std::unordered_map<int, Component> mapC;
  int id = 0;
  std::unordered_set<int> setf;
  // 1 - transform the vector of features in a set;
  for(std::vector<int>::iterator itv = annotD.anGroups[anG].begin(); itv != annotD.anGroups[anG].end(); itv++)
    setf.insert(*itv);
  // 2 - Compute components
  while(setf.size() != 0) {
    Component comp;
    std::unordered_set<int> extraf;
    std::unordered_set<int>::iterator it = setf.begin();
    extraf.insert(*it);
    getComponent(setm, extraf, annotD, comp);
    // drop features if comp created
    if(comp.feature.size() > 0) {
      for(std::unordered_set<int>::iterator itf = comp.feature.begin(); itf != comp.feature.end(); itf++)
	setf.erase(*itf);
      mapC[id] = comp;
      id++;
    }
  }
  return mapC;
}



class outputAn {
 public:
  std::vector<int> features;
  std::unordered_map<int, std::string> an1;
  std::unordered_map<int, double> mass1;
  std::unordered_map<int, double> score1;
  std::unordered_map<int, std::string> an2;
  std::unordered_map<int, double> mass2;
  std::unordered_map<int, double> score2;
  std::unordered_map<int, std::string> an3;
  std::unordered_map<int, double> mass3;
  std::unordered_map<int, double> score3;
  std::unordered_map<int, std::string> an4;
  std::unordered_map<int, double> mass4;
  std::unordered_map<int, double> score4;
  std::unordered_map<int, std::string> an5;
  std::unordered_map<int, double> mass5;
  std::unordered_map<int, double> score5;
};

outputAn createoutputAnno(annotDF& annotdf) {
  outputAn outAn;
  outAn.features = annotdf.features;
  for(std::vector<int>::iterator itv = outAn.features.begin(); itv != outAn.features.end(); itv++) {
    outAn.an1[*itv] = "NA";
    outAn.an2[*itv] = "NA";
    outAn.an3[*itv] = "NA";
    outAn.an4[*itv] = "NA";
    outAn.an5[*itv] = "NA";
    outAn.mass1[*itv] = 0;
    outAn.mass2[*itv] = 0;
    outAn.mass3[*itv] = 0;
    outAn.mass4[*itv] = 0;
    outAn.mass5[*itv] = 0;
    outAn.score1[*itv] = 0;
    outAn.score2[*itv] = 0;
    outAn.score3[*itv] = 0;
    outAn.score4[*itv] = 0;
    outAn.score5[*itv] = 0;
  }
  return outAn;
}



outputAn getAnnotation(Rcpp::DataFrame dfclique, Rcpp::DataFrame dfaddlist, int topmassf = 1, int topmasstotal = 10,
		       unsigned int sizeanG = 20, double tol = 0.00001, double filter = 0.0001, double emptyS = -6,
		       double newmass = -10.0, bool normalizeScore = true) {
  std::vector<int> topAn;
  // 1 - read ordered data frame of features and masses from R
  annotDF annotdf = readDF(dfclique);
  // 2 - read data frame with adduct list and adduct information from R
  rawadList rList = readrawList(dfaddlist);
  std::vector<double> vScore = getScoreAddlist(rList);
  //for(std::vector<double>::reverse_iterator itv = vScore.rbegin(); itv != vScore.rend(); itv++)
  // 3 - obtain all adducts and mass candidates for the data frame
  annotData annotD = getannotData(rList, annotdf, tol, filter );
  // 4 - create an object for putting the results of annotation
  outputAn outAn = createoutputAnno(annotdf);
  // 5 - find annotation for all annotation groups in this clique
  for(std::unordered_map<int, std::vector<int>>::iterator itg = annotD.anGroups.begin(); itg != annotD.anGroups.end(); itg++) {
    std::unordered_set<double> setm = getTopScoringMasses(annotD, itg->first, rList, topmassf, topmasstotal, emptyS);
    if(itg->second.size() > sizeanG) { // in case that there are a lot of features in this annotation group, separate components
      std::unordered_map<int, Component> components = getSeparateComp(setm, annotD, itg->first);
      // get all annotations
      for(std::unordered_map<int, Component>::iterator itc = components.begin(); itc != components.end(); itc++) {
	std::unordered_map<int, Annotation> annotations = AnnotateMassLoop(annotD, components[itc->first].feature, rList, 
									   components[itc->first].mass, emptyS);
	// get the id of the top five annotations
	topAn = sortAnnotations(annotations, 5);
	// 6 - Now put this annotations in the output object
	if(0 < topAn.size()) {
	  // normalize the scores
	  if(normalizeScore == true) {
	    for( std::vector<int>::iterator itv = topAn.begin(); itv!= topAn.end(); itv++)
	      normalizeAnnotation(annotations[*itv], vScore, newmass);
	  }
	  // annotation 1
	  for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita1 = annotations[topAn[0]].annotation.begin();
	      ita1 != annotations[topAn[0]].annotation.end(); ita1++) {
	    outAn.score1[ita1->first] = annotations[topAn[0]].score;
	    outAn.an1[ita1->first] = ita1->second.second;
	    outAn.mass1[ita1->first] = ita1->second.first;
	  }
	}
	if(1 < topAn.size()) {
	  // annotation 2
	  for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita2 = annotations[topAn[1]].annotation.begin();
	      ita2 != annotations[topAn[1]].annotation.end(); ita2++) {
	    outAn.score2[ita2->first] = annotations[topAn[1]].score;
	    outAn.an2[ita2->first] = ita2->second.second;
	    outAn.mass2[ita2->first] = ita2->second.first;
	  }
	}
	if(2 < topAn.size()) {
	  // annotation 3
	  for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita3 = annotations[topAn[2]].annotation.begin();
	      ita3 != annotations[topAn[2]].annotation.end(); ita3++) {
	    outAn.score3[ita3->first] = annotations[topAn[2]].score;
	    outAn.an3[ita3->first] = ita3->second.second;
	    outAn.mass3[ita3->first] = ita3->second.first;
	  }
	}
	if(3 < topAn.size()) {
	  // annotation 4
	  for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita4 = annotations[topAn[3]].annotation.begin();
	      ita4 != annotations[topAn[3]].annotation.end(); ita4++) {
	    outAn.score4[ita4->first] = annotations[topAn[3]].score;
	    outAn.an4[ita4->first] = ita4->second.second;
	    outAn.mass4[ita4->first] = ita4->second.first;
	  }
	}
	if(4 < topAn.size()) {
	  // annotation 5
	  for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita5 = annotations[topAn[4]].annotation.begin();
	      ita5 != annotations[topAn[4]].annotation.end(); ita5++) {
	    outAn.score5[ita5->first] = annotations[topAn[4]].score;
	    outAn.an5[ita5->first] = ita5->second.second;
	    outAn.mass5[ita5->first] = ita5->second.first;
	  }
	}
      }
    } else { // not necessary to separate in components, annotate all features in the anGroup
      std::unordered_set<int> setf;
      for(std::vector<int>::iterator itgf = annotD.anGroups[itg->first].begin();
	  itgf != annotD.anGroups[itg->first].end(); itgf++)
	setf.insert(*itgf);
      std::unordered_map<int, Annotation> annotations = AnnotateMassLoop(annotD, setf , rList, setm, emptyS);
      for(std::unordered_map<int, Annotation>::iterator itex1 = annotations.begin(); itex1 != annotations.end(); itex1++) {
	for(std::unordered_map<int, std::pair<double, std::string>>::iterator itex2 = annotations[itex1->first].annotation.begin();
	    itex2 != annotations[itex1->first].annotation.end(); itex2++) {
	}
      }
      topAn = sortAnnotations(annotations, 5);
      if(0 < topAn.size()) {
	if(normalizeScore == true) {
	  for( std::vector<int>::iterator itv = topAn.begin(); itv!= topAn.end(); itv++)
	    normalizeAnnotation(annotations[*itv], vScore, newmass);
	}
	// annotation 1
	for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita1 = annotations[topAn[0]].annotation.begin();
	    ita1 != annotations[topAn[0]].annotation.end(); ita1++) {
	  outAn.score1[ita1->first] = annotations[topAn[0]].score;
	  outAn.an1[ita1->first] = ita1->second.second;
	  outAn.mass1[ita1->first] = ita1->second.first;
	}
      }
      if(1 < topAn.size()) {
	// annotation 2
	for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita2 = annotations[topAn[1]].annotation.begin();
	    ita2 != annotations[topAn[1]].annotation.end(); ita2++) {
	  outAn.score2[ita2->first] = annotations[topAn[1]].score;
	  outAn.an2[ita2->first] = ita2->second.second;
	  outAn.mass2[ita2->first] = ita2->second.first;
	}
      }
      if(2 < topAn.size()) {
	// annotation 3
	for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita3 = annotations[topAn[2]].annotation.begin();
	    ita3 != annotations[topAn[2]].annotation.end(); ita3++) {
	  outAn.score3[ita3->first] = annotations[topAn[2]].score;
	  outAn.an3[ita3->first] = ita3->second.second;
	  outAn.mass3[ita3->first] = ita3->second.first;
	}
      }
      if(3 < topAn.size()) {
	// annotation 4
	for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita4 = annotations[topAn[3]].annotation.begin();
	    ita4 != annotations[topAn[3]].annotation.end(); ita4++) {
	  outAn.score4[ita4->first] = annotations[topAn[3]].score;
	  outAn.an4[ita4->first] = ita4->second.second;
	  outAn.mass4[ita4->first] = ita4->second.first;
	}
      }
      if(4 < topAn.size()) {
	// annotation 5
	for(std::unordered_map<int, std::pair<double, std::string>>::iterator ita5 = annotations[topAn[4]].annotation.begin();
	    ita5 != annotations[topAn[4]].annotation.end(); ita5++) {
	  outAn.score5[ita5->first] = annotations[topAn[4]].score;
	  outAn.an5[ita5->first] = ita5->second.second;
	  outAn.mass5[ita5->first] = ita5->second.first;
	}
      }
      setf.clear();
    }
  }
  return outAn;
}

#endif
