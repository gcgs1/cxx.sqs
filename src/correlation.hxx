/*
 correlation.hxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <boost/numeric/ublas/vector.hpp>
#include <vector>

#include "neighbour.hxx"
#include "supercell.hxx"

namespace ublas = boost::numeric::ublas;

class Correlation {

public:

  Correlation();

  void Calculate(const Neighbour& neighbour_table,
                 Supercell& supercell);

  double ErrorFunction(double x);

  long number;
  std::vector<double> pair_clusters;
  std::vector<long> pair_count;
  std::vector<double> pair_correlations;
  ublas::vector<double> errors;

};
