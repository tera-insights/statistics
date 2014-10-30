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

// Interface for implementations of distribution that can perform ML estimation
// of its parameters given sufficient statistics.

// Additional support is provided for mixtures models of said distribution.

#if !defined _CLUS_DISTRIBUTION_H_
#define _CLUS_DISTRIBUTION_H_

#include <string>
#include <iostream>
#include <math.h>
#include <armadillo>

#include "regressor.h"

namespace CLUS {
class Distribution {
 protected:
  // The dimension of the distribution being represented.
  int dim;

  // The weight of this distribution in the mixture model.
  double weight;

 public:
  // Simplest constructor possible.
  Distribution(int dim, double weight = 1)
      : dim(dim),
        weight(weight) {
  }

  Distribution(const Distribution& other)
      : dim(other.dim),
        weight(other.weight) {
  }

  // Creates an object that can do the regression for this distribution if given the input
  virtual Regressor* CreateRegressor() = 0;

  // Returns true if weight is 0, in which case the distribution should not be used.
  inline bool HasZeroWeight() {
    return weight == 0;
  }

  // Returns the normalized input.
  virtual arma::vec Normalize(const arma::vec& x) {
    return arma::vec(x);
  }

  // Returns the probability that a point lies in this distribution.
  virtual double PDF(const arma::vec& x) = 0;

  // Updates the sufficient statistics given a new point and probability.
  virtual void Update(const arma::vec& x, double prob) = 0;

  // Estimates the parameters from the sufficient statistics using ML.
  virtual double Estimate() = 0;

  // Random initialization of the distribution
  virtual void RandomDistribution(int NrClusters) = 0;
};
}
#endif // _CLUS_DISTRIBUTION_H
