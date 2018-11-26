/*
 main.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <functional>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/numeric.hpp>

namespace ublas = boost::numeric::ublas;
using namespace std;

class Neighbour {
  
public:
  
  Neighbour(unsigned long nbas,
	    ublas::matrix<double> plat,
	    ublas::matrix<double> basis);
  
  void GetSiteList(double rmax);
  
private:

  unsigned long NumberOfAtoms;
  ublas::matrix<double> LatticeVectors;
  ublas::matrix<double> BasisVectors;
  ublas::matrix<double> Sites;

};

Neighbour::Neighbour(unsigned long nbas,
		     ublas::matrix<double> plat,
		     ublas::matrix<double> basis) {
  NumberOfAtoms = nbas;
  LatticeVectors = plat;
  BasisVectors = basis;
}

void Neighbour::GetSiteList(double rmax) {
  ublas::vector<long> rcut(3);
  ublas::vector<long> sext(3);
  for (long i = 0; i < 3; i++) {
      ublas::matrix_row<ublas::matrix<double>> p(LatticeVectors,i);
      rcut(i) = static_cast<long>(rmax/ublas::norm_2(p)+2);
      sext(i) = 2*rcut(i) + 1;
  }
  long product = boost::accumulate(sext,1,std::multiplies<long>());
  long smax = NumberOfAtoms*product;
  ublas::matrix<double> sites(smax,3);
  ublas::matrix_row<ublas::matrix<double>> p0(LatticeVectors,0);
  ublas::matrix_row<ublas::matrix<double>> p1(LatticeVectors,1);
  ublas::matrix_row<ublas::matrix<double>> p2(LatticeVectors,2);
  long icnt = 0;
  for (long i = 0; i < NumberOfAtoms; i++) {
    ublas::matrix_row<ublas::matrix<double>> b(BasisVectors,i);
    for (long j = -rcut(0); j <= rcut(0); j++) {
          for (long k = -rcut(1); k <= rcut(1); k++) {
              for (long l = -rcut(2); l <= rcut(2); l++) {
                  for (int a = 0; a < 3; a++) {
                    sites(icnt,a) = b(a) + j*p0(a) + k*p1(a) + l*p2(a);
                  }
                  icnt++;
		  cout << icnt << endl;
	      }
	  }
    }
  }
  Sites = sites;
  cout << Sites << endl;
}


int main(void) {

  unsigned long nbas = 1;
  ublas::matrix<double> plat(3,3);
  ublas::matrix<double> basis(nbas,3);
  double rmax = 2.0;

  /* Set Basis Vectors */
  basis(0,0) = 0.0; basis(0,1) = 0.0; basis(0,2) = 0.0;
  
  /* Set Lattice Vectors */
  plat(0,0) =  0.5; plat(0,1) =  0.5; plat(0,2) = -0.5;
  plat(1,0) = -0.5; plat(1,1) =  0.5; plat(1,2) =  0.5;
  plat(2,0) =  0.5; plat(2,1) = -0.5; plat(2,2) =  0.5;

  Neighbour Nhbr = Neighbour(nbas,plat,basis);
  Nhbr.GetSiteList(rmax);
    
  return 0;
}
