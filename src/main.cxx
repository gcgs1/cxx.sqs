/*
 main.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>

#include "clusters.hxx"
#include "neighbour.hxx"
#include "utils.hxx"

namespace ublas = boost::numeric::ublas;

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

  for (auto i = 0; i < clusters.Number; i++) {
    std::cout << clusters.PairClusters[i] << ' ';
    std::cout << clusters.PairCount[i] << '\n';
  }
      
  return 0;
}
