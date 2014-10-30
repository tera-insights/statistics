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

#ifndef _CLUS_BINARYOBLIQUESPLITTER_H
#define _CLUS_BINARYOBLIQUESPLITTER_H

#include <cmath>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <jsoncpp/json/json.h>

#include "general.h"
#include "binarysplitter.h"
#include "splitpointcomputation.h"
#include "statisticsgatherers.h"

namespace CLUS {
// Splitter that is able to produce oblique splits, i.e. splits of the form
// a1 * X1 + ... > a0
class BinaryObliqueSplitter: public BinarySplitter {
 protected:
  // number of datapoints used to produce the split
  int N;

  // statistics for discrete variables
  std::vector<ProbabilisticBinomialStatistics>  discreteStatistics;

  // statistics for continuous variables
  MultidimNormalStatistics continuousStatistics;

  // is the partitioning pure? (class labels completely separated)
  bool purePart;

public:
  // Default constructor
  BinaryObliqueSplitter()
      : BinarySplitter(),
        N(0),
        purePart(false) {
  }

  /** Construct object when dimensionality is known
    @param DDomainSize   vector of domain sizes for discrete variables
    @param CsplitDim     number of split continuous variables
    @param RegDim      number of regressor attributes
  */
  BinaryObliqueSplitter(const arma::ivec DDomainSize, int CsplitDim, int RegDim)
      : BinarySplitter(DDomainSize, CsplitDim, RegDim),
        N(0),
        purePart(false) {
  }

  // Copy constructor
  // @param aux   object to be copied
  BinaryObliqueSplitter(const BinaryObliqueSplitter& aux)
      : BinarySplitter(aux),
        N(aux.N),
        discreteStatistics(aux.discreteStatistics),
        continuousStatistics(aux.continuousStatistics),
        purePart(aux.purePart) {
  }

  // which of the three methods is used for computation of separating hyperplane
  enum MultidimNormalStatistics::SeparationType CSepHypType;

  // Should be called after ComputeSplitVariable().
  //
  // @return true if the children given by branch will be splitted in the future
  bool MoreSplits(int branch, int Min_no_datapoints) {
    if (purePart)
      return false;
    else if (branch == 0)
      return continuousStatistics.GetS_P1() >= Min_no_datapoints;
    else
      return continuousStatistics.GetS_P2() >= Min_no_datapoints;
  }

  /** Compute the probability to take the left branch
    @param Dvars    discrete inputs
    @param Cvars    continuous inputs
  */
  double ProbabilityLeft(const arma::ivec Dvars, const arma::vec Cvars) {
    if (ProbabilityLeftPrivate(Dvars, Cvars) > 0.5)
      return 1.0;
    else
      return 0.0;
  }

  int ChooseBranch( const arma::ivec Dvars, const arma::vec Cvars) {
    if (ProbabilityLeftPrivate(Dvars, Cvars) > 0.5)
      return 0;
    else
      return 1;
  }

  void InitializeSplitStatistics(void) {
<<<<<<< HEAD
    std::cout << "Initializing split stats. csplitDim: " << csplitDim << " regDim: " << regDim << std::endl;
=======
>>>>>>> 871f53673e4a2ae8c56097105abf22e59e27a05a
    discreteStatistics.resize(dsplitDim);
    for (int i = 0; i < dsplitDim; i++)
      discreteStatistics[i].ResetDomainSize(dDomainSize[i]);
    continuousStatistics.Resize(csplitDim+regDim);
  }

  void UpdateSplitStatistics(const arma::ivec Dvars, const arma::vec Cvars,
                             double p1I, double p2I, double probability) {
    N++;

    double p1 = p1I * probability;
    double p2 = p2I * probability;

    // update the discrete var statistics
    for (int i = 0; i < dsplitDim; i++)
      discreteStatistics[i].UpdateStatisticsP(Dvars[i], p1, p2);

    continuousStatistics.UpdateStatistics(Cvars, p1, p2);
  }

  void Merge(const BinaryObliqueSplitter& other) {
    for (int i = 0; i < dsplitDim; i++)
      discreteStatistics[i].Merge(other.discreteStatistics[i]);
    continuousStatistics.Merge(other.continuousStatistics);
  }

  void DeleteTemporaryStatistics(void) {
    // free the space ocupied by the suficient statistics
    discreteStatistics.resize(0);
    continuousStatistics.Resize(0);
  }

  int ComputeSplitVariable(int type = MultidimNormalStatistics::LDA) {
    double maxgini;

    if (N == 0)
      goto error;

    // compute the best oblique split
    maxgini = continuousStatistics.ComputeGiniGain(type);
    SplitVariable = continuousStatistics.GetSplitVariable();

    // look for a cathegorical split that is better

    for (int i = 0; i < dsplitDim; i++) {
      double gini = discreteStatistics[i].ComputeGiniGain();

      //cout << "Var: " << i << " gini=" << gini << endl;

      if (gini > maxgini) {
        maxgini = gini;
        SplitVariable = i;
      }
    }

    if (maxgini == 0)
      goto error; // nobody makes any improvement

    if (SplitVariable >= 0)
      splitSetProbability = discreteStatistics[SplitVariable].GetProbabilitySet();
    else
      SeparatingHyperplane = continuousStatistics.GetSeparatingHyperplane();

    // anything is hopefully done at this point

    // if we are very close to maximum achievable gini we do no more splits in future
    if (fabs(maxgini - continuousStatistics.MaxGini()) < TNNearlyZero)
      purePart=true;

    return 0;

    error:
    std::cout << "Error encountered, killing node" << std::endl;

    return -1;
  }

 private:
  /** Computes the probability to follow the left branch
    @param Dvars    discrete inputs
    @param Cvars    continuous inputs
  */
  double ProbabilityLeftPrivate(const arma::ivec Dvars, const arma::vec Cvars) {
    // use the fact that SeparatingHyperplane is normalized
    if (SplitVariable == -1) {
      // split on a hyperplane
      double crit = SeparatingHyperplane[0]
                  + dot(Cvars,
                        SeparatingHyperplane.subvec(1, csplitDim + regDim));
      return 1.0 - PValueNormalDistribution(0.0, 1.0, crit);
    } else  if (SplitVariable <= -2) {
      // split on a continuous variable
      double crit = SeparatingHyperplane[0]
                  + Cvars[-SplitVariable - 2]
                  * SeparatingHyperplane[-SplitVariable - 1];
      return  1.0 - PValueNormalDistribution(0.0, 1.0, crit);
    } else {
      // split on a discrete variable
      int value = Dvars[SplitVariable];
      return splitSetProbability[value];
    }
  }
 public:
  Json::Value ToJson() const {
    Json::Value ret;
    if (SplitVariable == -1) {
      ret["type"] = "oblique";
      for (auto value : SeparatingHyperplane)
        ret["plane"].append(value);
    } else if (SplitVariable <= -2) {
      ret["type"] = "continuous";
      ret["variable"] = -SplitVariable - 2;
      ret["line"].append(SeparatingHyperplane[0]);
      ret["line"].append(SeparatingHyperplane[-SplitVariable - 1]);
    } else {
      ret["type"] = "discrete";
      ret["variable"] = SplitVariable;
      for (auto value : splitSetProbability)
        ret["set"].append(value);
    }
    return ret;
  }
};
}

#endif // _CLUS_BINARYOBLIQUESPLITTER_H
