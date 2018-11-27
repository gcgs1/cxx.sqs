/*
 main.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/io.hpp>

#include "neighbour.hxx"
#include "utils.hxx"

int main(void) {

  unsigned long nbas = 1;
  //unsigned long nbas = 2;
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
  
  Neighbour Nhbr = Neighbour(nbas,plat,basis);
  Nhbr.SetInverseLattice(iplat,xbasis);
  Nhbr.GetNhbrList(rmax);

  std::cout << Nhbr.Total << std::endl;
  std::cout << Nhbr.Pointers << std::endl;
  
  return 0;
}
