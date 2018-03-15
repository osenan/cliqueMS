// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef NETWORKClIQUEMSR_H
#define NETWORKCLIQUEMSR_H

#include <RcppArmadilloExtensions/sample.h>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cmath>


// One include file from Boost
#include <boost/functional/hash.hpp>

typedef std::pair <int, int> Keyd; // key consiting of edge1 and edge 2
typedef boost::hash<std::pair<int, int> > Hashd;
typedef std::unordered_map<Keyd, double, Hashd> Edges;
typedef Edges::const_iterator edgeKey; // iterator for edges
typedef std::unordered_map<Keyd, bool, Hashd> Edgeclique; // unorderred map category to indicate if the edge is inside a clique

// New Function createEdges but starting from a Data Frame
Edges createEdges (Rcpp::DataFrame netdf)
{
  Edges edges; // dictionary of edges
  double weight;
  int node1;
  int node2;
  Rcpp::NumericVector vnode1 = netdf["node1"];
  Rcpp::NumericVector vnode2 = netdf["node2"];
  Rcpp::NumericVector vweigth = netdf["weight"];
  for(unsigned int index = 0; index < vnode1.size(); index++) {
    node1 = vnode1[index];
    node2 = vnode2[index];
    weight = vweigth[index];
    edges[std::make_pair(node1, node2)] = weight;
  }
  return edges;
}
  
// Function createNodes, with node x you can acces clique value y
std::unordered_map<int, int> createNodes (Edges edges)
{
  std::unordered_map<int, int> nodes; // declare nodes unordered_map
  std::unordered_map<int, int>::const_iterator got1; // iterator for node1 for nodes object
  std::unordered_map<int, int>::const_iterator got2; // iterator for node2 for nodes object
  int node1;
  int node2;
  for ( edgeKey it(edges.begin()); it != edges.end(); it++ ) // iterate over edge keys
    { node1 = it->first.first; 
      node2 = it->first.second;
      got1 = nodes.find(node1); // check if node1 is in keys of node map
      got2 = nodes.find(node2); // check if node2 is in keys of node map
      if( got1 == nodes.end() ) // if is not in key, add it and initialize clique attribute
	nodes[node1] = node1; 
      if( got2 == nodes.end() ) // if is not in key, add it and initialize clique attribute
	nodes[node2] = node2; 
    }
  return nodes;
}

// Function createNeighbours for creating vector of neighbors accesed by key node
std::unordered_map<int, std::vector<int> > createNeighbors (Edges edges)
{
  std::unordered_map<int, std::vector<int> > neighbors; // initialize neighbours
  int node1;
  int node2;
  for ( edgeKey it(edges.begin()); it != edges.end(); it++ ) // iterate over edge keys
    { node1 = it->first.first; 
      node2 = it->first.second;
      neighbors[node1].push_back(node2);
      neighbors[node2].push_back(node1);
    }
  return neighbors;
}

// Function createEdgeclique for creating boolean variable to indicate if this edge is inside or not a clique
Edgeclique createEdgecliques (Edges edges)
{
  Edgeclique edgecliques; // dictionary of edges containing if they belong or not to a clique
  int node1;
  int node2;
  bool entry =false; // if edges are inside cliques is True, that is relevant for log likelihood calculation
  for ( edgeKey it(edges.begin()); it != edges.end(); it++ ) // iterate over edges pair values
    {
      node1 = it->first.first;
      node2 = it->first.second;
      edgecliques[std::make_pair(node1, node2)] = entry; // at the beggining, all edges are outside of cliques
    }
  return edgecliques;
}

// Function createCliques to create a unordered map of key clique name with a vector of the nodes of that clique
std::unordered_map<int, std::vector<int> > createCliques (std::unordered_map<int, int> nodes)
{
  std::unordered_map<int, std::vector<int> > cliques; // initialize cliques
  int node;
  int cliqueV;
  for (std::unordered_map<int, int>::iterator it = nodes.begin() ; it != nodes.end(); ++it)
    {
      node = it->first;
      cliqueV = it->second;
      cliques[node].push_back(cliqueV);
    }
  return cliques;
}

class Network {
public:
  Edges edges;
  std::unordered_map<int, int> nodes;
  std::unordered_map<int, std::vector<int> > neighbors;
  std::unordered_map<int, std::vector<int> > cliques;
  Edgeclique edgecliques;
  Edges logedges; // log value of weigth powered to some exponent
  Edges minuslogedges; // 1 - log value of weight powered to some exponent
};

// Function to initialize network
Network createNetwork(Rcpp::DataFrame netdf, double exp = 2)
{
  Network net;
  net.edges = createEdges(netdf);
  net.nodes = createNodes(net.edges);
  net.neighbors = createNeighbors(net.edges);
  net.cliques = createCliques(net.nodes);
  net.edgecliques = createEdgecliques(net.edges);
  int node1;
  int node2;
  double weight, logpower, minuslogpower;
  for ( edgeKey it(net.edges.begin()); it != net.edges.end(); it++ )
    {
      node1 = it->first.first;
      node2 = it->first.second;
      weight = net.edges[std::make_pair(node1,node2)];
      logpower = log10(pow(weight,exp));
      minuslogpower = log10(1 -pow(weight,exp));
      net.logedges[std::make_pair(node1,node2)] = logpower;
      net.minuslogedges[std::make_pair(node1,node2)] = minuslogpower;
    }
  return net;
}

// Function to calculate log-likelihood from a network
double logltotal(Network& net)
{
  double inside = 0.0, outside = 0.0, logl = 0.0;
  int node1;
  int node2;
  std::pair < int, int> edgepair;
  for ( Edges::const_iterator it(net.edges.begin()); it != net.edges.end(); it++ ) // iterate over edges
    {
      node1 = it->first.first;
      node2 = it->first.second;
      edgepair = std::make_pair(node1,node2);
      if( net.edgecliques[edgepair] == true) {
	  inside += net.logedges[edgepair];
	  //	  	  std::cout << "inside " << inside << "\n";
      } else {
	outside += net.minuslogedges[edgepair];
	//std::cout << "outside"  << outside << " \n";
      }
    }
  logl = inside + outside;

  return logl;
}

// Nodelogl class to return an object with the change in logl, the new edges inside clique and the oldeges than will become outside edges
class Nodelogl { 
public:
  double logl;
  int newnode;
  std::vector<std::pair<int,int> > newedges;
  std::vector<std::pair<int,int> > oldedges;
};

void sortEdge(std::pair<int, int>& edge)
{
  int ifirst = edge.first;
  int isecond = edge.second;
  if( ifirst > isecond ) { // if the key to edges value is in the wrong order, correct the order: always s1 < s2
    edge.first = isecond;
    edge.second = ifirst;
  }
}

// Function to compute the change of likelihood if we add node2 to clique of node1
Nodelogl calcNodelogl(Network& net, int node1, int node2)
{
  Nodelogl nResult;
  double logl = -1.0, logl_change = 0.0, logl_before = 0.0;
  bool complete = true;
  double newlinks_change = 0.0, nolinks_change = 0.0, newlinks_before = 0, nolinks_before = 0;
  std::vector<std::pair<int,int> > newEdges;
  std::vector<std::pair<int,int> > oldEdges;
  int clique1 = net.nodes[node1];
  int clique2 = net.nodes[node2];
  std::vector<int> nclique1 = net.cliques[clique1];
  nclique1.erase(std::find(nclique1.begin(),nclique1.end(),node1)); // remove the node1 apart from the other nodes of its clique
  std::vector<int> nclique2 = net.cliques[clique2];
  nclique2.erase(std::find(nclique2.begin(),nclique2.end(),node2)); // remove the node2 apart from the other nodes of its clique
  if(nclique1.size() > 0) {
    for(std::vector<int>::iterator it1 = nclique1.begin(); it1 != nclique1.end(); it1++) {
      std::pair<int,int> edge1 = std::make_pair(*it1,node2);
      sortEdge(edge1);
      if(net.logedges.find(edge1) == net.logedges.end()) {
	complete = false;
	break;
      } else {
	newlinks_change += net.logedges[edge1]; // edges that now will be part of the clique
	newlinks_before += net.minuslogedges[edge1]; // this edges were before outside cliques
	newEdges.push_back(edge1);
      }
    }
  }
  if(nclique2.size() > 0) {
    for(std::vector<int>::iterator it2 = nclique2.begin(); it2 != nclique2.end(); it2++) {
      std::pair<int,int> edge2 = std::make_pair(*it2,node2);
      sortEdge(edge2);
      nolinks_change += net.minuslogedges[edge2]; // this edges now will be outside cliques
      nolinks_before += net.logedges[edge2]; // this edges between node2 and its old clique members were inside cliques before
      oldEdges.push_back(edge2);
    }
  }
  if(complete == true) {
    std::pair<int,int> edge = std::make_pair(node1,node2);
    sortEdge(edge);
    newlinks_change += net.logedges[edge];
    newlinks_before += net.minuslogedges[edge];
    newEdges.push_back(edge);
    logl_change = newlinks_change + nolinks_change;
    logl_before = newlinks_before + nolinks_before;
    logl = logl_change - logl_before;
  }
  nResult.logl = logl;
  nResult.newnode = node1; // initialize results values in case there is a break
  nResult.newedges = newEdges;
  nResult.oldedges = oldEdges;
  return nResult;
}

// function to move the current node to another clique
double reassignNode(Network& net, int node, double logl)
{
  Nodelogl maxchange;
  maxchange.logl = 0;
  if(net.neighbors[node].size() > 0) { // reassign this node if it has neighbors
    int ownclique = net.nodes[node];
    std::unordered_set<int> diffcliques;
    for(std::vector<int>::iterator itn = net.neighbors[node].begin(); itn!= net.neighbors[node].end(); itn++) {
      int cliquecandidate = net.nodes[*itn];
      if(cliquecandidate != ownclique) 
	if(diffcliques.find(cliquecandidate) == diffcliques.end() )
	  diffcliques.insert(cliquecandidate);
    }
    if(diffcliques.size() > 0) {
      for(std::unordered_set<int>::iterator itc = diffcliques.begin(); itc!= diffcliques.end(); itc++) {
	int node2 = net.cliques[*itc][0]; // first node in the clique candidate
	Nodelogl nodeChange = calcNodelogl(net, node2, node); // logl change if we move node to clique of node2
	if(nodeChange.logl > maxchange.logl) // we search for the max change in logl bigger than 0
	  maxchange = nodeChange;
      }
    }
    if(maxchange.logl > 0) { // if there is a positive change in logl by moving a node, now execute this change
      for(std::vector<std::pair<int,int> >::iterator itnew = maxchange.newedges.begin();
	  itnew != maxchange.newedges.end(); itnew++)
	net.edgecliques[*itnew] = true; // change makesnew edges are inside cliques
      for(std::vector<std::pair<int,int> >::iterator itold = maxchange.oldedges.begin();
	  itold != maxchange.oldedges.end(); itold++)
	net.edgecliques[*itold] = false; // change puts edges outside cliques
      int newclique = net.nodes[maxchange.newnode];
      net.nodes[node] = newclique; // now the clique of node is the clique of node2
      net.cliques[ownclique].erase(std::find(net.cliques[ownclique].begin(),
					     net.cliques[ownclique].end(),node)); // remove the node from clique nodes of old clique
      net.cliques[newclique].push_back(node); // add node to the list of nodes of the new clique
    }
  }
  double newlogl = logl + maxchange.logl;
  return newlogl;
}    

Nodelogl calcCliquelogl(Network& net, int clique1, int clique2)
{
  Nodelogl cResult;
  double logl = -1.0, logl_change = 0.0, logl_before = 0.0;
  bool complete = true;
  std::vector<std::pair<int,int> > newEdges;
  std::vector<std::pair<int,int> > oldEdges;
  for(std::vector<int>::iterator it1 = net.cliques[clique1].begin(); it1 != net.cliques[clique1].end(); it1++) {
    for(std::vector<int>::iterator it2 = net.cliques[clique2].begin(); it2 != net.cliques[clique2].end(); it2++) {
      std::pair<int,int> edge1 = std::make_pair(*it1,*it2);
      sortEdge(edge1);
      if(net.logedges.find(edge1) == net.logedges.end()) {
	complete = false;
	goto endloglCalc; // exit the two loops if one link between the two cliques does not exist
      } else {
	logl_change += net.logedges[edge1]; // edges that now will be part of the clique
	logl_before += net.minuslogedges[edge1]; // this edges were before outside cliques
	newEdges.push_back(edge1);
      }
    }
  }
 endloglCalc: 
  if(complete == true)
    logl = logl_change -logl_before;
  cResult.newnode = clique1; // initialize results values in case there is a break
  cResult.logl = logl;
  cResult.newedges = newEdges;
  return cResult;
}

double meanClique(Network&net, int clique1, int clique2)
{
  double meanV = 0.0; double meanR = -1;
  double size = 0.0;
  bool complete = true;
  for(std::vector<int>::iterator it1 = net.cliques[clique1].begin(); it1 != net.cliques[clique1].end(); it1++) {
    for(std::vector<int>::iterator it2 = net.cliques[clique2].begin(); it2 != net.cliques[clique2].end(); it2++) {
      std::pair<int,int> edge1 = std::make_pair(*it1,*it2);
      sortEdge(edge1);
      if(net.logedges.find(edge1) == net.logedges.end()) {
	complete = false;
	goto endMeancalc; // exit the two loops if one link between the two cliques does not exist
      } else {
	double weight = net.edges[edge1];
	meanV += pow(weight,2.0);
	size++;
      }
    }
  }
 endMeancalc:
  if( complete == true)
    meanR = meanV/size;
  return meanR;
}


double reassignClique(Network& net, int clique, double logl) {
  double loglchange = 0; // the change in logl if we accept joining clique to another clique
  int node = net.cliques[clique][0];
  std::unordered_set<int> diffcliques;
  for(std::vector<int>::iterator itn = net.neighbors[node].begin(); itn!= net.neighbors[node].end(); itn++) {
    int cliquecandidate = net.nodes[*itn];
      if(cliquecandidate != clique) 
	if(diffcliques.find(cliquecandidate) == diffcliques.end() )
	  diffcliques.insert(cliquecandidate);
  }
  if(diffcliques.size() >0 ) {
    double maxmean = 0.0;
    int maxclique;
    for(std::unordered_set<int>::iterator itc = diffcliques.begin(); itc!= diffcliques.end(); itc++) {
      double meancandidate = meanClique(net, clique, *itc);
      if(meancandidate > maxmean) {
	maxclique = *itc;
	maxmean = meancandidate;
      }
    }
    Nodelogl loglcandidate = calcCliquelogl(net, clique, maxclique);
    if(loglcandidate.logl > 0) {
      loglchange = loglcandidate.logl;
      for(std::vector<std::pair<int,int> >::iterator itedge = loglcandidate.newedges.begin();
	  itedge != loglcandidate.newedges.end(); itedge++)
	net.edgecliques[*itedge] = true; // change makes new edges be inside cliques
      for(std::vector<int>::iterator itnode = net.cliques[clique].begin();
	  itnode != net.cliques[clique].end(); itnode++) {
	net.nodes[*itnode] = maxclique; // change clique value of old clique nodes to the new
	net.cliques[maxclique].push_back(*itnode); // add nodes of old clique to the list of nodes of the new clique
      }
      net.cliques.erase(clique); // delete old clique because it is empty
    }
  }
  double loglReturn = logl + loglchange;
  return loglReturn;
}

Rcpp::NumericVector csample_integer(Rcpp::NumericVector x, int size, bool replace,
				    Rcpp::NumericVector prob = Rcpp::NumericVector::create()) {
  Rcpp::NumericVector ret = Rcpp::RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

std::vector<double> itReassign(Network& net, double tol, double logl) {
  double currentlogl = logl;
  std::vector<double> loglResult;
  loglResult.push_back(currentlogl);
  Rcpp::NumericVector allnodes;
  Rcpp::NumericVector randallnodes;
  for(std::unordered_map<int,int>::iterator n = net.nodes.begin(); n != net.nodes.end(); n++) 
    allnodes.push_back(n->first); // insert allnodes values in vector of nodes
  // the order of nodes has to be random
  randallnodes = csample_integer(allnodes, allnodes.size(), false);
  for(Rcpp::NumericVector::iterator itv = randallnodes.begin(); itv != randallnodes.end(); itv++) {
    currentlogl = reassignNode(net, *itv, currentlogl);
    loglResult.push_back(currentlogl);
  }
  double firstlogl = loglResult[0];
  double diff = 1  - abs(currentlogl/firstlogl); // diference in log likelihood after one complete round of node reassignments
  int rcount = 1; // counter of the number of rounds
  while( diff > tol ) {
    double firstlogl = loglResult.back();
    randallnodes = csample_integer(allnodes, allnodes.size(), false);
    for(Rcpp::NumericVector::iterator itv = randallnodes.begin(); itv != randallnodes.end(); itv++) {
      currentlogl = reassignNode(net, *itv, currentlogl); // move nodes to different cliques
      loglResult.push_back(currentlogl); // store results of change in logl
    }
    diff = 1  - abs(currentlogl/firstlogl); // diference in log likelihood after one complete round of node reassignments
    rcount++;
  }
  Rcpp::Rcout << "Kernighan-Lin done with " << rcount << " rounds\n";
  return loglResult;
}


std::vector<double> aggregateANDkernighan(Network& net, double tol, int step) {
  double currentlogl = logltotal(net);
  std::vector<double> loglResult;
  loglResult.push_back(currentlogl);
  // round 1
  Rcpp::NumericVector allnodes;
  Rcpp::NumericVector randallnodes;
  for(std::unordered_map<int,int>::iterator n = net.nodes.begin(); n != net.nodes.end(); n++)
      allnodes.push_back(n->first); // insert allnodes values in vector of nodes
  randallnodes = csample_integer(allnodes, allnodes.size(), false);
  int scount = 1; // counter of the number of rounds that are clique joining, it starts with 1
  int tcount = 1; // total number of rounds
  for(Rcpp::NumericVector::iterator itv = randallnodes.begin(); itv != randallnodes.end(); itv++) {
    int cliquec = net.nodes[*itv]; // clique that will be joined to another clique
    currentlogl = reassignClique(net, cliquec, currentlogl);
    loglResult.push_back(currentlogl);
    scount++;
    tcount++;
    if(scount == step) {
      currentlogl = reassignNode(net, *itv, currentlogl);
      loglResult.push_back(currentlogl);
      scount = 1;
      tcount++;
    }
  }
  double firstlogl = loglResult[0];
  double diff = 1  - abs(currentlogl/firstlogl); // diference in log likelihood after one complete round of node reassignments
    // rest of rounds
  while( diff > tol ) {
    firstlogl = loglResult.back();
    randallnodes = csample_integer(allnodes, allnodes.size(), false);
    for(Rcpp::NumericVector::iterator itw = randallnodes.begin(); itw != randallnodes.end(); itw++) {
      int cliquecw = net.nodes[*itw]; // clique that will be joined to another clique
      currentlogl = reassignClique(net, cliquecw, currentlogl);
      loglResult.push_back(currentlogl);
      scount++;
      tcount++;
      if(scount == step) {
	currentlogl = reassignNode(net, *itw, currentlogl);
	loglResult.push_back(currentlogl);
	scount = 1;
	tcount++;
      }
    }
    diff = 1  - abs(currentlogl/firstlogl);
  }
  Rcpp::Rcout <<"Aggregate cliques done, with " << tcount << " rounds\n";
  // Kernighan-Lin after aggregation of cliques
  std::vector<double> loglLast = itReassign(net, tol, currentlogl);
  for(std::vector<double>::iterator itf = loglLast.begin() ; itf != loglLast.end(); itf++)
    loglResult.push_back(*itf);
  return loglResult;
}

#endif
