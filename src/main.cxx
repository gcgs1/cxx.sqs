/*
 main.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <functional>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/numeric.hpp>

#include "utils.hxx"

class Neighbour {
  
public:
  
  Neighbour(unsigned long nbas,
	    ublas::matrix<double> plat,
	    ublas::matrix<double> basis);

  void SetInverseLatticeVectors(ublas::matrix<double> iplat);
  void GetSiteList(double rmax);
  void GetNhbrList(double rmax);
  ublas::vector<long> LinkedListBlocks(double rmax);
  
private:

  unsigned long NumberOfAtoms;
  ublas::matrix<double> LatticeVectors;
  ublas::matrix<double> BasisVectors;
  ublas::matrix<double> Sites;
  ublas::vector<long> SupercellExtension;
  ublas::matrix<double> InverseLatticeVectors;
  
};

Neighbour::Neighbour(unsigned long nbas,
		     ublas::matrix<double> plat,
		     ublas::matrix<double> basis) {
  NumberOfAtoms = nbas;
  LatticeVectors = plat;
  BasisVectors = basis;
}

void Neighbour::SetInverseLatticeVectors(ublas::matrix<double> iplat) {
  InverseLatticeVectors = iplat;
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
	}
      }
    }
  }
  Sites = sites;
  SupercellExtension = sext;
}

ublas::vector<long> Neighbour::LinkedListBlocks(double rmax) {

  /* Decompose Lattice Vector Matrix into Lattice Vectors */
  ublas::matrix_row<ublas::matrix<double>> p1(LatticeVectors,0);
  ublas::matrix_row<ublas::matrix<double>> p2(LatticeVectors,1);
  ublas::matrix_row<ublas::matrix<double>> p3(LatticeVectors,2);

  /* Cross Product of Lattice Vectors */
  ublas::vector<double> c1 = cross3(p2,p3);
  ublas::vector<double> c2 = cross3(p3,p1);
  ublas::vector<double> c3 = cross3(p1,p2);

  ublas::vector<double> w(3);
  w(0) = ublas::inner_prod(p1,c1)/ublas::norm_2(c1);
  w(1) = ublas::inner_prod(p2,c2)/ublas::norm_2(c2);
  w(2) = ublas::inner_prod(p3,c3)/ublas::norm_2(c3);

  ublas::vector<long> blks(3);
  for (long i = 0; i < 3; i++) {
    blks(i) = static_cast<long>(SupercellExtension(i)*w(i)/rmax);
    if (blks(i) < 3) blks(i) = 3;
  }

  return blks;
}

void Neighbour::GetNhbrList(double rmax) {

  /* Calculate Site List */
  Neighbour::GetSiteList(rmax);

  /* Linked-List Block Size */
  ublas::vector<long> blks = Neighbour::LinkedListBlocks(rmax);

  ublas::vector<long> nhbr_tot(650*NumberOfAtoms);
  
  unsigned long nsum = 1;
  for (long int i = 0; i < NumberOfAtoms; i++) {
    nhbr_tot(i) = 0;
    ublas::matrix_row<ublas::matrix<double>> b(BasisVectors,i);
    ublas::vector<double> bx = ublas::prod(b,InverseLatticeVectors);
    std::cout << bx << std::endl;
    //bx = matmul(s_cell%basis(i,:),s_cell%iplat)
  }
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

  ublas::matrix<double> iplat = inv3(plat);
  
  Neighbour Nhbr = Neighbour(nbas,plat,basis);
  Nhbr.SetInverseLatticeVectors(iplat);
  Nhbr.GetNhbrList(rmax);
    
  return 0;
}
