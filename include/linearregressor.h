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

// Implements a linear regressor, for which y = x^T * beta + k. Note that the
// intercept is stored as its own field in order to avoid padding input.

#if !defined _CLUS_LINEARREGRESSOR_H
#define _CLUS_LINEARREGRESSOR_H

#include <armadillo>
#include <iostream>
#include <jsoncpp/json/json.h>

#include "regressor.h"

namespace CLUS {
class LinearRegressor: public Regressor {
 protected:
  arma::vec beta;

  double intercept;

 public:
  // Default constructor.
  LinearRegressor()
      : beta(),
        intercept() {
  }

  // Simple constructor, in which the intercept is the first term of beta.
  LinearRegressor(arma::vec& beta)
      : beta(beta.subvec(1, beta.n_elem - 1)),
        intercept(beta(0)) {
  }

  // Explicit constructor.
  LinearRegressor(arma::vec& beta, double intercept)
      : beta(beta),
        intercept(intercept) {
  }

  // Copy Constructor.
  LinearRegressor(const LinearRegressor& other)
      : beta(other.beta),
        intercept(other.intercept) {
  }

  virtual ~LinearRegressor() {}

  // Performs the prediction by taking the inner product of x with beta.
  virtual inline double Predict(const arma::vec& x) const override {
    return dot(x, beta) + intercept;
  }

  virtual inline Regressor* clone() const override {
    return new LinearRegressor(*this);
  }

  virtual Json::Value ToJson() const override {
    Json::Value ret;
    ret.append(intercept);
    for (auto value : beta)
      ret.append(value);
    return ret;
  }
};
}

#endif // _CLUS_LINEARREGRESSOR_H
