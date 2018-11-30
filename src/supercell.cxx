/*
 supercell.cxx

 Copyright (c) 2018 Guy Skinner
 
 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "supercell.hxx"
#include "utils.hxx"

namespace ublas = boost::numeric::ublas;

Supercell::Supercell(ublas::vector<long> sext,
                     long primitive_number_of_atoms,
                     ublas::matrix<double> primitive_lattice_vectors,
                     ublas::matrix<double> primitive_basis_vectors,
		     double concentration) {

  supercell_extension = sext;
  primitive_number_of_atoms_ = primitive_number_of_atoms;
  primitive_lattice_vectors_ = primitive_lattice_vectors;
  primitive_basis_vectors_ = primitive_basis_vectors;
  concentration_ = concentration;

  number_of_atoms = product(supercell_extension)*primitive_number_of_atoms_;

  lattice_vectors = primitive_lattice_vectors_;
  for (auto i = 0; i < 3; i++) {
    long s = supercell_extension(i);
    for (auto j = 0; j < 3; j++) {
      double pij = primitive_lattice_vectors_(i,j);
      lattice_vectors(i,j) = s*pij;
    }
  }

  ublas::matrix<double> basis(number_of_atoms,3);
  auto iat = 0;
  for (auto i = 0; i < primitive_number_of_atoms_; i++) {
    for (auto j = 0; j < supercell_extension(0); j++) {
      for (auto k = 0; k < supercell_extension(1); k++) {
        for (auto l = 0; l < supercell_extension(2); l++) {
          for (auto a = 0; a < 3; a++) {
	    basis(iat,a) = primitive_basis_vectors_(i,a)
	                 + j*primitive_lattice_vectors_(0,a)
                         + k*primitive_lattice_vectors_(1,a)
	                 + l*primitive_lattice_vectors_(2,a);
          }
          ublas::matrix_row<ublas::matrix<double>> b(basis_vectors,iat);
          iat++;
        }
      }
    }
  }
  
  auto solute = concentration_*number_of_atoms;
  ublas::vector<long> tmp(number_of_atoms);
  for (auto i = 0; i < number_of_atoms; i++) {
    if (i < solute) {
      tmp(i) = -1;
    } else {
      tmp(i) = 1;
    }
  }

  basis_vectors = basis;
  pointers = tmp;
  
}
