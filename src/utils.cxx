/*
 utils.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <cmath>
#include "utils.hxx"

double det3(ublas::matrix<double> a) {
  /* Fast Determinant of 3x3 Matrix */
  double det = a(0,0)*(a(1,1)*a(2,2) - a(1,2)*a(2,1))
             - a(0,1)*(a(1,0)*a(2,2) - a(1,2)*a(2,0))
             + a(0,2)*(a(1,0)*a(2,1) - a(1,1)*a(2,0));
  return det;
}

ublas::vector<double> cross3(ublas::vector<double> a,
                             ublas::vector<double> b) {

  /* Cross Product for Three-Dimensional Vectors */
  if (a.size() != 3 or b.size() != 3) {
    std::cout << "Exit(1): Vector dimensions not compatible for cross product\n";
    exit(1);
  }
  ublas::vector<double> c(3);
  c(0) = a(1)*b(2) - a(2)*b(1);
  c(1) = a(2)*b(0) - a(0)*b(2);
  c(2) = a(0)*b(1) - a(1)*b(0);
  return c;
}

ublas::matrix<double> inv3(ublas::matrix<double> a) {

  /* Matrix Inversion for Three-Dimensional Vectors */
  ublas::matrix<double> ia(a.size1(),a.size2());
  double det = det3(a);
  
  ia(0,0) = a(1,1)*a(2,2) - a(1,2)*a(2,1);
  ia(0,1) = a(0,2)*a(2,1) - a(0,1)*a(2,2);
  ia(0,2) = a(0,1)*a(1,2) - a(0,2)*a(1,1);
  
  ia(1,0) = a(1,2)*a(2,0) - a(1,0)*a(2,2);
  ia(1,1) = a(0,0)*a(2,2) - a(0,2)*a(2,0);
  ia(1,2) = a(0,2)*a(1,0) - a(0,0)*a(1,2);
  
  ia(2,0) = a(1,0)*a(2,1) - a(1,1)*a(2,0);
  ia(2,1) = a(0,1)*a(2,0) - a(0,0)*a(2,1);
  ia(2,2) = a(0,0)*a(1,1) - a(1,0)*a(0,1);
  
  return ia/det;
}

long linear_search(const std::vector<double>& vec,
		   double a,
		   double tol) {

  for (long i = 0; i < vec.size(); i++) {
    double value = fabs(vec[i] - a);
    if (value < tol) {
      return i;
    }
  }
  
  return vec.size();
}
