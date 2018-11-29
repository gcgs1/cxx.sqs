/*
 clusters.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <numeric>
#include <vector>

#include "clusters.hxx"
#include "utils.hxx"

namespace ublas = boost::numeric::ublas;

Clusters::Clusters() {
}

void Clusters::FindClusters(const Neighbour& neighbour_table,
                            ublas::matrix<double> basis) {

  /* Loop Over Neighbours to Find, Number and Sort Pairs */
  ublas::matrix<double> table = neighbour_table.Sites;
  std::vector<double> pairs = {};
  std::vector<long> count = {};

  /* Assume Equivalent Atoms */
  ublas::matrix_row<ublas::matrix<double>> b(basis,0);

  for (auto j = 0; j < neighbour_table.Total(0); j++) {
    auto ptr = neighbour_table.Pointers(j);
    ublas::matrix_row<ublas::matrix<double>> site(table,ptr);
    double r = ublas::norm_2(site - b);
    long index = linear_search(pairs,r,1e-6);
    if (index == pairs.size()) {
      pairs.insert(pairs.begin(),r);
      count.insert(count.begin(),1);
    } else {
      count[index]++;
    }
  }

  /* Sorting */
  std::vector<std::size_t> indices(pairs.size());
  std::iota(indices.begin(),indices.end(),0);
  std::sort(indices.begin(), indices.end(),
            [&pairs](std::size_t left, std::size_t right) {
              return pairs[left] < pairs[right];
            });

  std::vector<double> sort_pairs = pairs;
  std::vector<long> sort_count = count;

  Number = pairs.size();
  for (auto i = 0; i < Number; i++) {
    sort_pairs[i] = pairs[indices[i]];
    sort_count[i] = count[indices[i]];
  }

  PairClusters = sort_pairs;
  PairCount = sort_count;

}