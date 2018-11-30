/*
 neighbours.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "neighbour.hxx"
#include "utils.hxx"

Neighbour::Neighbour(unsigned long nbas,
                     ublas::matrix<double> plat,
                     ublas::matrix<double> basis) {
  NumberOfAtoms = nbas;
  LatticeVectors = plat;
  BasisVectors = basis;
}

Neighbour::~Neighbour() {
}

void Neighbour::SetInverseLattice(ublas::matrix<double> iplat,
                                  ublas::matrix<double> xbasis) {
  InverseLatticeVectors = iplat;
  InverseBasisVectors = xbasis;
}

void Neighbour::GetSiteList(double rmax) {

  /* Create Site List */
  ublas::vector<long> rcut(3);
  ublas::vector<long> sext(3);
  
  for (auto i = 0; i < 3; i++) {
    ublas::matrix_row<ublas::matrix<double>> p(LatticeVectors,i);
    rcut(i) = static_cast<long>(rmax/ublas::norm_2(p)+2);
    sext(i) = 2*rcut(i) + 1;
  }
  
  SupercellExtension = sext;
  /* Linked-List Block Size */
  Neighbour::LinkedListBlocks(rmax);
  
  ublas::vector<double> cent(3);
  for (auto i = 0; i < 3; i++) cent(i) = 0.5*(sext(i) - 1);
  Centre = cent;
  
  //long product = boost::accumulate(sext,1,std::multiplies<long>());
  long smax = NumberOfAtoms*product(sext);
  long blks = boost::accumulate(Blocks,1,std::multiplies<long>());
  
  ublas::matrix<double> sites(smax,3);
  ublas::vector<long> basptr(smax);
  ublas::vector<long> linkl(smax);
  ublas::vector<long> head(blks);
  ublas::vector<double> sitex(3);

  /* Make Sure Zero Vectors */
  head.clear();
  linkl.clear();
  
  ublas::matrix_row<ublas::matrix<double>> p0(LatticeVectors,0);
  ublas::matrix_row<ublas::matrix<double>> p1(LatticeVectors,1);
  ublas::matrix_row<ublas::matrix<double>> p2(LatticeVectors,2);
  
  long icnt = 0;
  for (auto i = 0; i < NumberOfAtoms; i++) {
    ublas::matrix_row<ublas::matrix<double>> b(BasisVectors,i);
    ublas::matrix_row<ublas::matrix<double>> xb(InverseBasisVectors,i);
    for (auto j = -rcut(0); j <= rcut(0); j++) {
      for (auto k = -rcut(1); k <= rcut(1); k++) {
        for (auto l = -rcut(2); l <= rcut(2); l++) {
          for (auto a = 0; a < 3; a++) {
            sites(icnt,a) = b(a) + j*p0(a) + k*p1(a) + l*p2(a);
            sitex(a) = xb(a);
          }
          basptr(icnt) = i;
          sitex(0) += static_cast<double>(j);
          sitex(1) += static_cast<double>(k);
          sitex(2) += static_cast<double>(l);
          unsigned long icell = FindCell(sitex);
          linkl(icnt) = head(icell);
          head(icell) = icnt;
          icnt++;
        }
      }
    }
  }

  Sites = sites;
  BasisPointers = basptr;
  Head = head;
  LinkedList = linkl;
  
}

void Neighbour::LinkedListBlocks(double rmax) {

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
  for (auto i = 0; i < 3; i++) {
    blks(i) = static_cast<long>(SupercellExtension(i)*w(i)/rmax);
    if (blks(i) < 3) blks(i) = 3;
  }

  Blocks = blks;

}

void Neighbour::GetNhbrList(double rmax) {

  /* Calculate Site List */
  Neighbour::GetSiteList(rmax);

  ublas::vector<long> tot(NumberOfAtoms);
  ublas::vector<long> ptrs(650*NumberOfAtoms);
  
  unsigned long nsum = 0;
  for (auto i = 0; i < NumberOfAtoms; i++) {
    tot(i) = 0;
    ublas::matrix_row<ublas::matrix<double>> b(BasisVectors,i);
    ublas::matrix_row<ublas::matrix<double>> bx(InverseBasisVectors,i);
    unsigned long icell = FindCell(bx);
    for (auto imx = -1; imx <= 1; imx++) {
      for (auto imy = -Blocks(0); imy <= Blocks(0); imy += Blocks(0)) {
        for (auto imz = -Blocks(0)*Blocks(1);
             imz <= Blocks(0)*Blocks(1);
             imz += Blocks(0)*Blocks(1)) {
          long jcell = icell + imx + imy + imz;
          long j = Head(jcell);
          while (j != 0) {
            ublas::matrix_row<ublas::matrix<double>> site(Sites,j);
            ublas::vector<double> ra = site - b;
            double r = ublas::norm_2(ra);
            if (r <= rmax and r > 0.0) {
              ptrs(nsum) = j;
              tot(i)++;
              nsum++;
            }
            j = LinkedList(j);
          }
        }
      }
    }
  }

  Total = tot;
  Pointers = ptrs;

}

unsigned long Neighbour::FindCell(ublas::vector<double> sitex){

  /* Find Block that Atom Belongs In */
  ublas::vector<long> ixs(3);
  ublas::vector<double> scent = sitex + Centre;

  for (auto i = 0; i < 3; i++) {
    ixs(i) = static_cast<long>(Blocks(i)*scent(i)/SupercellExtension(i));
    if (ixs(i) > Blocks(i)) ixs(i) = Blocks(i)-1;
  }

  unsigned long icell = ixs(0) + Blocks(0)*ixs(1) + Blocks(1)*Blocks(0)*ixs(2);
  return icell;

}
