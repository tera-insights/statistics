/*
Copyright (c) 2014, Tera Insights, LLC.
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

#ifndef _CLUS_REGRESSION_TREE_H_
#define _CLUS_REGRESSION_TREE_H_

#include <armadillo>
#include <cstddef>
#include <limits>

// Interface class for regression trees.

class RegressionTree {
public:
  using SizeType = std::size_t;

  using ContinuousValue = double;
  using ContinuousVector = arma::vec;
  using DiscreteValue = arma::sword;
  using DiscreteVector = arma::ivec;

  virtual ~RegressionTree() {}

	virtual void StartLearningRound(void) = 0;
	virtual void LearnSample(
      const DiscreteVector& dVars,
      const ContinuousVector& cVars,
      const double probability = 1.0) = 0;
	virtual bool StopLearningRound(void) = 0;

	virtual void StartPruningRound(void) = 0;
	virtual void PruningSample(
      const DiscreteVector& dVars,
      const ContinuousVector& cVars,
      const ContinuousValue y,
      const double probability = 1.0) = 0;
  virtual bool StopPruningRound(void) = 0;

	virtual ContinuousValue Infer(
      const DiscreteVector& dVars,
      const ContinuousVector& cVars,
      const double infer_threshold,
      const SizeType maxNodeID = std::numeric_limits<SizeType>::max()) = 0;
};

#endif // _CLUS_REGRESSION_TREE_H_
