// [[Rcpp::depends(BH)]]
#include "networkCliqueMSR.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::DataFrame returnCliques(Rcpp::DataFrame netdf, double tol = 0.000001) {
  Network net = createNetwork(netdf);
  double logl = logltotal(net);
  std::vector<int> vnode;
  std::vector<int> vclique;
  Rcpp::Rcout << "Beggining value of logl is " << logl << " \n";
  std::vector<double> loglList = aggregateANDkernighan(net, tol, logl);
  for(std::unordered_map<int,int>::iterator n = net.nodes.begin(); n != net.nodes.end(); n++) {
    vnode.push_back(n->first);
    vclique.push_back(n->second);
  }
  double loglfinal = logltotal(net);
  Rcpp::Rcout << "Finishing value of logl is " << loglfinal << " \n";
  return Rcpp::DataFrame::create(Rcpp::Named("node") = vnode, Rcpp::Named("clique") = vclique);
}
