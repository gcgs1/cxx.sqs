/*
 utils.hxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

double det3(ublas::matrix<double> a);

ublas::vector<double> cross3(ublas::vector<double> a,
                             ublas::vector<double> b);

ublas::matrix<double> inv3(ublas::matrix<double> a);

long linear_search(const std::vector<double>& vec,
		   double a,
		   double tol);
