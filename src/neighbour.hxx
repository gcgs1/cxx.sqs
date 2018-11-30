/*
 neighbours.hxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

class Neighbour {
  
public:
  
  Neighbour(unsigned long nbas,
            ublas::matrix<double> plat,
            ublas::matrix<double> basis);

  ~Neighbour();
  
  void SetInverseLattice(ublas::matrix<double> iplat,
                         ublas::matrix<double> xbasis);
  void GetSiteList(double rmax);
  void GetNhbrList(double rmax);
  void LinkedListBlocks(double rmax);
  unsigned long FindCell(ublas::vector<double> sitex);
  
  ublas::vector<long> Total;
  ublas::vector<long> Pointers;
  ublas::vector<long> BasisPointers;
  ublas::matrix<double> Sites;
  
private:
  
  long NumberOfAtoms;
  ublas::matrix<double> LatticeVectors;
  ublas::matrix<double> BasisVectors;
  ublas::vector<long> SupercellExtension;
  ublas::matrix<double> InverseLatticeVectors;
  ublas::matrix<double> InverseBasisVectors;
  ublas::vector<long> Blocks;
  ublas::vector<double> Centre;
  ublas::vector<long> Head;
  ublas::vector<long> LinkedList;
  
};
