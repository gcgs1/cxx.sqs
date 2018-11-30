/*
 supercell.hxx

 Copyright (c) 2018 Guy Skinner
 
 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

class Supercell {

public:

  Supercell(ublas::vector<long> supercell_extension,
            long primitive_number_of_atoms,
            ublas::matrix<double> primitive_lattice_vectors,
            ublas::matrix<double> primitive_basis_vectors,
	    double concentration);

  ublas::vector<long> supercell_extension;
  long number_of_atoms;
  ublas::matrix<double> lattice_vectors;
  ublas::matrix<double> basis_vectors;
  ublas::vector<long> pointers;
  
private:

  long primitive_number_of_atoms_;
  ublas::matrix<double> primitive_lattice_vectors_;
  ublas::matrix<double> primitive_basis_vectors_;
  double concentration_;

};
