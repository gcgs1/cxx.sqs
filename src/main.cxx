/*
 main.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "neighbour.hxx"
#include "utils.hxx"

namespace ublas = boost::numeric::ublas;

class Clusters {

public:

  Clusters();
  void FindClusters(const Neighbour& neighbour_table,
		    ublas::matrix<double> basis);

  std::vector<double> PairClusters;
  std::vector<long> PairCount;
  
private:

};

Clusters::Clusters() {
}

void Clusters::FindClusters(const Neighbour& neighbour_table,
			    ublas::matrix<double> basis) {

  /* Loop Over Neighbours to Find, Number and Sort Pairs */
  ublas::matrix<double> table = neighbour_table.Sites;
  std::vector<double> pairs = {};
  std::vector<long> count = {};
  
  /* Assume Equivalent Atoms */
  ublas::matrix_row<ublas::matrix<double>> b(basis,0);
    
  for (auto j = 0; j < neighbour_table.Total(0); j++) {
    auto ptr = neighbour_table.Pointers(j);
    ublas::matrix_row<ublas::matrix<double>> site(table,ptr);
    double r = ublas::norm_2(site - b);
    long index = linear_search(pairs,r,1e-6);
    if (index == pairs.size()) {
      pairs.insert(pairs.begin(),r);
      count.insert(count.begin(),1);
    } else {
      count[index]++;
    }
  }

  /* Sorting */
  std::vector<std::size_t> indices(pairs.size());
  std::iota(indices.begin(),indices.end(),0);
  std::sort(indices.begin(), indices.end(),
	    [&pairs](std::size_t left, std::size_t right) {
	      return pairs[left] < pairs[right];
	    });

  std::vector<double> sort_pairs = pairs;
  std::vector<long> sort_count = count;
  
  for (auto i = 0; i < pairs.size(); i++) {
    sort_pairs[i] = pairs[indices[i]];
    sort_count[i] = count[indices[i]];
  }

  PairClusters = sort_pairs;
  PairCount = sort_count;
  
}

int main(void) {

  //unsigned long nbas = 1;
  unsigned long nbas = 1;
  ublas::matrix<double> plat(3,3);
  ublas::matrix<double> basis(nbas,3);
  double rmax = 2.0;

  /* Set Basis Vectors */
  basis(0,0) = 0.0; basis(0,1) = 0.0; basis(0,2) = 0.0;

  /*
  double q = 1.62;
  basis(0,0) = 0.0; basis(0,1) = 0.0; basis(0,2) = 0.0;
  basis(1,0) = 0.5; basis(1,1) = std::sqrt(3)/6; basis(1,2) = 0.5*q;
  */

  std::cout << basis << std::endl;
  
  /* Set Lattice Vectors */
  plat(0,0) =  0.5; plat(0,1) =  0.5; plat(0,2) = -0.5;
  plat(1,0) = -0.5; plat(1,1) =  0.5; plat(1,2) =  0.5;
  plat(2,0) =  0.5; plat(2,1) = -0.5; plat(2,2) =  0.5;

  /*
  plat(0,0) =  1.0; plat(0,1) =              0.0; plat(0,2) = 0.0;
  plat(1,0) = -0.5; plat(1,1) = 0.5*std::sqrt(3); plat(1,2) = 0.0;
  plat(2,0) =  1.0; plat(2,1) =              0.0; plat(2,2) =   q;
  */

  std::cout << plat << std::endl;

  //exit(0);
  
  ublas::matrix<double> iplat = inv3(plat);
  ublas::matrix<double> xbasis = ublas::prod(basis,iplat);
  
  Neighbour neighbour_table = Neighbour(nbas,plat,basis);
  neighbour_table.SetInverseLattice(iplat,xbasis);
  neighbour_table.GetNhbrList(rmax);
  
  Clusters clusters = Clusters();
  clusters.FindClusters(neighbour_table,basis);
    
  return 0;
}
