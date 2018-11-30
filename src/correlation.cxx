/*
 correlation.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <numeric>
#include <vector>

#include "correlation.hxx"
#include "neighbour.hxx"
#include "supercell.hxx"
#include "utils.hxx"

Correlation::Correlation() {
}

void Correlation::Calculate(const Neighbour& neighbour_table,
                            Supercell& supercell) {

  ublas::matrix<double> table = neighbour_table.Sites;

  std::vector<double> pairs = {};
  std::vector<long> count = {};
  std::vector<double> corrs = {};

  auto nsum = 0;

  for (auto i = 0; i < supercell.number_of_atoms; i++) {
    ublas::matrix_row<ublas::matrix<double>> ib(supercell.basis_vectors,i);
    auto ip = supercell.pointers(i);
    for (auto j = 0; j < neighbour_table.Total(i); j++) {
      auto ptr = neighbour_table.Pointers(nsum+j);
      auto jbsp = neighbour_table.BasisPointers(ptr);
      auto jp = supercell.pointers(jbsp);

      ublas::matrix_row<ublas::matrix<double>> jb(table,ptr);

      double r = ublas::norm_2(jb-ib);
      double corr = static_cast<double>(ip*jp);
      long index = linear_search(pairs,r,1e-6);

      if (index == pairs.size()) {
        pairs.insert(pairs.begin(),r);
        count.insert(count.begin(),1);
        corrs.insert(corrs.begin(),corr);
      } else {
        count[index]++;
        corrs[index] += corr;
      }
    }
    nsum += neighbour_table.Total(i);
  }

  /* Normalize */
  for (auto i = 0; i < pairs.size(); i++) {
    corrs[i] /= count[i];
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
  std::vector<double> sort_corrs = corrs;

  number = pairs.size();
  for (auto i = 0; i < number; i++) {
    sort_pairs[i] = pairs[indices[i]];
    sort_count[i] = count[indices[i]];
    sort_corrs[i] = corrs[indices[i]];
  }

  pair_clusters = sort_pairs;
  pair_count = sort_count;
  pair_correlations = sort_corrs;

}

double Correlation::ErrorFunction(double x) {

  double error = (2*x-1)*(2*x-1);
  ublas::vector<double> efs(number);

  for (auto i = 0; i < number; i++) {
    efs(i) = std::fabs(error-pair_correlations[i]);
  }

  errors = efs;
  return ublas::sum(efs)/number;

}
