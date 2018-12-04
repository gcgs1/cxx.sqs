/*
 main.cxx

 Copyright (c) 2018 Guy Skinner

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <boost/format.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/range/algorithm.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

#include "correlation.hxx"
#include "neighbour.hxx"
#include "supercell.hxx"
#include "utils.hxx"

namespace ublas = boost::numeric::ublas;

int main(void) {

  unsigned long nbas = 2;
  ublas::matrix<double> plat(3,3);
  ublas::matrix<double> basis(nbas,3);
  double rmax = 2.0;

  /* Set Basis Vectors */
  double isqrt2 = 1/std::sqrt(2);
  double isqrt3 = 1/std::sqrt(3);
  double isqrt6 = 1/std::sqrt(6);

  basis(0,0) =      0.0; basis(0,1) =          0.0; basis(0,2) =            0.0;
  basis(1,0) =   isqrt6; basis(1,1) =       isqrt2; basis(1,2) =         isqrt3;

  plat(0,0)  = 2*isqrt6; plat(0,1)  =          0.0; plat(0,2)  = std::sqrt(3)/6;
  plat(1,0)  =      0.0; plat(1,1)  = std::sqrt(2); plat(1,2)  =            0.0;
  plat(2,0)  =      0.0; plat(2,1)  =          0.0; plat(2,2)  = std::sqrt(3)/2;

  /* Primitive I/O */
  std::cout << boost::format("Primitive Number of Atoms : %8d\n") % nbas;
  std::cout << boost::format("Pair Cut Off Radius       : %16.7f\n") % rmax;
  std::cout << boost::format("Primitive Lattice Vectors : %16.7f %16.7f %16.7f\n")
    % plat(0,0) % plat(0,1) % plat(0,2);
  std::cout << boost::format("                            %16.7f %16.7f %16.7f\n")
    % plat(1,0) % plat(1,1) % plat(1,2);
  std::cout << boost::format("                            %16.7f %16.7f %16.7f\n")
    % plat(2,0) % plat(2,1) % plat(2,2);

  for (auto i = 0; i < nbas; i++) {
    ublas::matrix_row<ublas::matrix<double>> b(basis,i);
    if (i == 0) {
      std::cout << boost::format("Primitive Basis Vectors   : %16.7f %16.7f %16.7f\n")
	% b(0) % b(1) % b(2);
    } else {
      std::cout << boost::format("                          : %16.7f %16.7f %16.7f\n")
        % b(0) % b(1) % b(2);
    }
  }
  
  /* Supercell Creation */
  ublas::vector<long> sext(3);
  sext(0) = 100; sext(1) = 60; sext(2) = 3;

  double concentration = 0.05;
  
  Supercell supercell = Supercell(sext,nbas,plat,basis,concentration);

  /* Supercell Inverse */
  ublas::matrix<double> isplat = inv3(supercell.lattice_vectors);
  ublas::matrix<double> xsbasis(supercell.number_of_atoms,3);
  for (auto i = 0; i < supercell.number_of_atoms; i++) {
    ublas::matrix_row<ublas::matrix<double>> b(supercell.basis_vectors,i);
    auto xb = prod(b,isplat);
    for (auto a = 0; a < 3; a++) {
      xsbasis(i,a) = xb(a);
    }
  }
  
  /* Supercell I/O */
  std::cout << boost::format("Supercell Number of Atoms : %8d\n")
    % supercell.number_of_atoms;
    
  std::cout << boost::format("Supercell Lattice Vectors : %16.7f %16.7f %16.7f\n")
    % supercell.lattice_vectors(0,0)
    % supercell.lattice_vectors(0,1)
    % supercell.lattice_vectors(0,2);
  std::cout << boost::format("                            %16.7f %16.7f %16.7f\n")
    % supercell.lattice_vectors(1,0)
    % supercell.lattice_vectors(1,1)
    % supercell.lattice_vectors(1,2);
  std::cout << boost::format("                            %16.7f %16.7f %16.7f\n")
    % supercell.lattice_vectors(2,0)
    % supercell.lattice_vectors(2,1)
    % supercell.lattice_vectors(2,2);
  
  /* Random Number Generator */
  boost::mt19937 mt(time(0));
  boost::uniform_int<> uni_dist;
  boost::variate_generator<boost::mt19937&,boost::uniform_int<>>
    generator(mt,uni_dist);

  long nitf = 100000;
  double term = 0.0001;

  std::cout << boost::format("Number of Iterations      : %8d\n") % nitf;
  std::cout << boost::format("Concentration             : %16.7f\n") % concentration;
  std::cout << boost::format("Perfect Correlation       : %16.7f\n")
    % ((2*concentration-1)*(2*concentration-1)); 
  std::cout << boost::format("SQS Accept Tolerance      : %16.7f\n\n") % term;

  Neighbour neighbour = Neighbour(supercell.number_of_atoms,
				  supercell.lattice_vectors,
				  supercell.basis_vectors);
  
  neighbour.SetInverseLattice(isplat,xsbasis);
  neighbour.GetNhbrList(rmax);

  double best = 1.0;
  long snum = 0;
  long bnum = 0;
  
  for (long n = 0; n < nitf; n++) { 
    
    boost::range::random_shuffle(supercell.pointers,generator);
        
    /* Correlation Function Calculation */
    Correlation correlation = Correlation();
    correlation.Calculate(neighbour,supercell);
    double result = correlation.ErrorFunction(concentration);

    bool bestq = (result < best);
    bool writeq = (result < term);
    if (bestq) {
      best = result;
      bnum = snum;
    }
    
    std::cout << boost::format("Iteration                 : %8d\n") % (n+1);
    std::cout << boost::format("Mean Error                : %16.7f\n")
      % result;
    std::cout << boost::format("Least Error      %8d : %16.7f\n") % bnum % best;
    if (writeq) {
      snum++;
      std::cout << "Accept\n";
    } else {
      std::cout << "Reject\n";
    }
    
    for (auto i = 0; i < correlation.number; i++) {
      std::cout << boost::format("  %16.7f  %5d  %16.7f  %16.7f\n")
	% correlation.pair_clusters[i]
	% (correlation.pair_count[i]/supercell.number_of_atoms)
	% correlation.pair_correlations[i]
	% correlation.errors[i];
    }
    std::cout << std::endl;

    /* File Writing */
    if (writeq) {
      std::ofstream outfile;
      boost::format fmt = boost::format("sqs.%d.out") % snum;
      outfile.open(fmt.str());
      outfile << boost::format("# special quasirandom structure\n");
      for (auto i = 0; i < supercell.number_of_atoms; i++) {
	long ptr = 1;
	if (supercell.pointers(i) == -1) ptr = 2;
	outfile << boost::format("%8d %22.15f %22.15f %22.15f\n")
	  % ptr
	  % supercell.basis_vectors(i,0)
	  % supercell.basis_vectors(i,1)
	  % supercell.basis_vectors(i,2);
      }
      outfile.close();
    }
    
  }
  
  return 0;

}
