/*

Copyright (c) 2003, Cornell University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   - Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.
   - Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
   - Neither the name of Cornell University nor the names of its
       contributors may be used to endorse or promote products derived from
       this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

*/

// This class implements a multi-dimensional normal distribution.
// For notational consistency, consult the following article:
// http://en.wikipedia.org/wiki/Multivariate_normal_distribution

#if !defined _CLUS_MULTINORMAL_H_
#define _CLUS_MULTINORMAL_H_

#include <math.h>
#include <armadillo>
#include <string>
#include <iostream>
#include <cmath>
#include <cassert>

#include "general.h"
#include "distribution.h"
#include "linearregressor.h"
#include "exceptions.h"

namespace CLUS {
class Multinormal : public Distribution {
 protected:
  // The mu parameter of the distribution. This can be estimated as the average
  // of the tuples.
  arma::vec mu;

  // The Cholesky decomposition of sigma, i.e. chol * chol^T = sigma. Because
  // all computations are performed using chol, sigma is not saved. Note that
  // the ML estimate of sigma is biased, which reflects on chol.
  arma::mat chol;

  // The coefficient of the pdf of the multinormal distribution, weighted
  // accordingly to account for the mixture model.
  double coef;

  // The number of tuples processed.
  int count;

  // The sum of the probability weights associated with each tuple.
  double sum_p;

  // The sum of tuples weighted by their probabilities.
  arma::vec sum_px;

  // Sum of the outer product of tuples with themselves weighted by probability.
  arma::mat sum_pxxT;

  // The radius of the distribution, used to generate nearby distributions.
  double radius;

  // The equation of the linear regressor, including a constant term.
  arma::vec line;

 public:
  // Basic Constructor.
  Multinormal(int dim = 0)
      : Distribution(dim),
        mu(dim + 1),
        chol(dim + 1, dim + 1),
        count(0),
        sum_p(0),
        sum_px(dim + 1, arma::fill::zeros),
        sum_pxxT(dim + 1, dim + 1, arma::fill::zeros),
        line(dim + 1) {
  }

  // Constructor of a multinormal distribution given its parameters
  Multinormal(int dim, double weight, arma::vec& mu, arma::mat& chol)
      : Distribution(dim, weight),
        mu(mu),
        chol(chol),
        count(0),
        sum_p(0),
        sum_px(dim + 1, arma::fill::zeros),
        sum_pxxT(dim + 1, dim + 1, arma::fill::zeros),
        line(dim + 1) {
    ComputeCoef();
  }

  // Copy constructor.
  Multinormal(const Multinormal& other)
      : Distribution(other),
        mu(other.mu),
        chol(other.chol),
        coef(other.coef),
        count(other.count),
        sum_p(other.sum_p),
        sum_px(other.sum_px),
        sum_pxxT(other.sum_pxxT),
        radius(other.radius),
        line(other.line) {
  }

 public:
  Regressor* CreateRegressor() {
    if (ComputeLine(chol))
      return new LinearRegressor(line);
    else
      return new LinearRegressor();
  }

  // Reverses the normalization transformation in the original space.
  void DenormalizeParameters(const Multinormal& parent) {
    mu = parent.chol * mu + parent.mu;
    chol = parent.chol * chol;
  }

  // Transforms a point belonging to this distribution such it appears to have
  // come from a standard multinormal distribution.
  inline arma::vec Normalize(const arma::vec& x) override {
    return chol.i() * (x - mu);
  }

  // Computes the PDF of the distribution for the given point.
  double PDF(const arma::vec& x) override {
    if (weight == 0)
      return 0;

    auto temp = chol.i() * (x - mu);
    double norm_sq = dot(temp, temp);
    return coef * std::exp(-norm_sq / 2.0);
  }

  // This is used to update the sufficient statistics of the distribution, which
  // will be used for the ML estimation of the parameters in UpdateParameters.
  void Update(const arma::vec& x, double prob) override {
    count++;
    sum_p += prob;
    sum_px += prob * x;
    sum_pxxT += prob * x * x.t();
  }

  void Merge(const Multinormal& other) {
    count += other.count;
    sum_p += other.sum_p;
    sum_px += other.sum_px;
    sum_pxxT += other.sum_pxxT;
  }

  // Computes chol and mu from the sufficient statistics and returns the mean
  // square of the change in mu from the previous estimate, with the initial
  // estimate of mu being treated as the origin.
  double Estimate() override {
    double distP = 0;
    int cholRet;
    arma::mat sigma;

    // if the cluster is dead do nothing
    if (weight == 0)
      return 0;

    if (sum_p <= (dim + 1)) {
      std::cerr << "Cluster was too small and was thrown out." << std::endl;
      weight = 0;
      distP = 1;
      goto cleanup;
    }

<<<<<<< HEAD
    // The parameters are estiamted
    mu = sum_px / sum_p;
    sigma = sum_pxxT / sum_p - mu * mu.t();
    distP = sum(square(mu - sum_px / sum_p));


    // The Cholesky decomposition of sigma is stored in chol. The transpose
    // is taken as Armadillo computes chol such that chol.t() * chol = sigma.
    cholRet = !arma::chol(chol, sigma);
    chol = chol.t();
    if (cholRet) {
      // If we got in here we have trouble
      // If the cluster does not depend on at least 2*(dim) points we just throw it out
      if (sum_p < 2 * (dim + 1)) {
        std::cerr << "Troublesome cluster too small to bother. We throw it out" << std::endl;
        weight = 0;
        distP = 1;
        goto cleanup;
      } else {
        std::cerr << "Matrix sigma is not positive definite:" << std::endl << sigma << std::endl;
        std::cerr << "chol:" << std::endl << chol;
        std::cerr << "weight=" << weight << std::endl;
        // fix the last value of chol
        if (cholRet == 2) {
          std::cout << "Fixing chol" << std::endl;
          chol(dim, dim) = 0; // only the last value is bad, we can fix it
          distP = -1; // this will make the returned value nan and terminate anything
        } else {
          // probably the cluster is a "squeezy"(has one less dimention in the input space).
          std::cerr << "The cluster is a squeezy, setting weight to 0." << std::endl;
          weight = 0;
          distP = 1 ;
          goto cleanup;
        }
      }
    }
    std::cout << "mu: " << mu.t();
    std::cout << "chol" << std::endl << chol << std::endl;

    ComputeCoef();

    cleanup:
    // Reset count, sum_p, sum_px, sum_pxxT
    ResetStatistics();

    return sqrt(distP) / (dim + 1); // how much the center has moved
  }

  void RandomDistribution(int NrClusters) override{
    const double smallness = 0.1;
    int d = dim + 1;
    const double sizeDist = 1 / pow(NrClusters, 1.0 / d);

    // reset weights to 1.
    weight = 1.0 / NrClusters;
    // choose mu random
    mu = mu.randu() * 2 - 1;

    do {
      // choose n randomly and build V, set of vectors orthogonal on n
      arma::vec n = arma::normalise(arma::randu<arma::vec>(d) * 2 - 1);

      // set Sigma as I-0.1n*n'
      // in this way Sigma has n as the eigenvector coresp to eightvalue 0.1
      // and d-1 other eigenvectors orth. on n with eigvals 1.0
      chol = sizeDist * arma::trimatu(arma::eye(d, d) -  smallness * n * n.t());
    } while(!arma::chol(chol, chol));
    chol = chol.t();

    ComputeCoef();
  }

  /** For hierarchies or trees we need to produce offsprings close to the parent
  @param NrClusters    total number of clusteres/distributions initialized
  @param parent      parent distribution in the hierachy
  */
  void RandomDistribution(int num_clusters, Multinormal& parent) {
    // reset weitht to 1.0. It might have been 0.
    weight = 1.0;

    // Create a distribution 1/num_clusters as small as the parent with a center
    // around the center of the parent
    radius = parent.radius / pow(num_clusters, 1.0 / dim);
    mu.subvec(0, dim - 1) = parent.mu + (2 * parent.radius - radius) * arma::randu<arma::vec>(dim);
    mu[dim] = parent.line[0] + dot(mu.subvec(0, dim - 1), parent.line.subvec(1, dim));

    // comment this line if the parent.C is used
    // mu[dim]=parent.mu[dim];

    // put sigma=I => chol=I*radius
    chol.fill(0);
    chol.diag().fill(radius);

    ComputeCoef();
  }

 protected:
  // Computes the coefficient of the PDF once weight and sigma have been determined.
  void ComputeCoef() {
    coef = weight * std::pow(2 * arma::datum::pi, -dim / 2) / det(chol);
  }

  /** Computes the equation of the line using Sigma=sum_pxxT.
  Must be called when sum_pxxT contains Sigma not only the partial sums for current implementation
  @param Sigma       the covariance matrix
  @param isCholesky  true if Sigma contains the Cholesky decomposition, false if original covariance matrix
  @return true if anything is fine, false if something is wrong and the computation is bogus.
  */
  bool ComputeLine(const arma::mat& Sigma) {
    std::cout << "Computing line using sigma:" << std::endl << Sigma << std::endl;
    if (weight == 0)
      return false;

    // use least square
    // form b, which is length dim + 1
    arma::vec b = Sigma.row(dim).t();
    b(dim) = 0.0;
    // compute in radius the trace of original S
    // suppose one of the eigenvalues is small.
    radius = sqrt(sum(square(Sigma.diag())) / dim);

    // solve for c = G \ b.
    // Note that last element is 0 because the last element of b is 0.
    arma::vec c = trimatu(Sigma).i() * b;

    // form line
    line[0] = mu[dim] - dot(mu, c);
    line.subvec(1, dim) = c.subvec(0, dim - 1);
    std::cout << "Result:" << std::endl << line << std::endl;
    return true;
  }

  void ResetStatistics() {
    count = 0;
    sum_p = 0;
    sum_px.zeros();
    sum_pxxT.zeros();
  }
};
}

#endif // _CLUS_MULTINORMAL_H_
