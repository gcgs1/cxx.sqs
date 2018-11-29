/*
 clusters.hxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "neighbour.hxx"

namespace ublas = boost::numeric::ublas;

class Clusters {

public:

  Clusters();
  void FindClusters(const Neighbour& neighbour_table,
                    ublas::matrix<double> basis);

  long Number;
  std::vector<double> PairClusters;
  std::vector<long> PairCount;

};
