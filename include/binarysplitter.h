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

#ifndef _CLUS_BINARYSPLITTER_H
#define _CLUS_BINARYSPLITTER_H

#include <iostream>
#include <armadillo>
#include <jsoncpp/json/json.h>

/** Base clases for all the splitters. Specifies the interface.

  Keeps track of the statistics for variables, decides what the
  split predicate should be and indicates split when input is
  provided.
*/
class BinarySplitter {
 protected:
  // the number of discrete split variables (nonregressers)
  int dsplitDim;

  // the number of continuous split variables (nonregressoers)
  int csplitDim;

  // the number of regressors
  int regDim;

  // list of discrete domain sizes
  const arma::ivec dDomainSize;

  /** Indicates on what variable this node splits on.
  If -1 we have an oblique split. */
  int SplitVariable;

  /** Contains the coeficients of the hyperplane that separates the 2 distributions
  in the input space. This is an oblique split */
  arma::vec SeparatingHyperplane;

  /** Contains the probability for a value to be in the left partition. In
  traditional classifiers is set to 1.0 or 0.0 */
  arma::vec splitSetProbability;

 public:
  // Default constructor
  BinarySplitter()
    : dDomainSize(),
      SeparatingHyperplane(),
      splitSetProbability() {
  }

  /** Construct object when dimensionality is known
  @param DDomainSize   vector of domain sizes for discrete variables
  @param CsplitDim     number of split continuous variables
  @param RegDim      number of regressor attributes
  */
  BinarySplitter(const arma::ivec DDomainSize, int CsplitDim, int RegDim)
    : dsplitDim(DDomainSize.n_elem),
      csplitDim(CsplitDim),
      regDim(RegDim),
      dDomainSize(DDomainSize),
      SeparatingHyperplane(),
      splitSetProbability() {
  }

  // Copy constructor
  // @param aux   object to be copied
  BinarySplitter(const BinarySplitter& aux)
    : dsplitDim(aux.dsplitDim),
      csplitDim(aux.csplitDim),
      regDim(aux.regDim),
      dDomainSize(aux.dDomainSize),
      SplitVariable(aux.SplitVariable),
      SeparatingHyperplane(aux.SeparatingHyperplane),
      splitSetProbability(aux.splitSetProbability) {
  }

  // Get the number of regressor variables
  inline int GetRegDim() {
    return regDim;
  }

  // Get the number of split variables
  inline int GetCSplitDim() {
    return csplitDim;
  }

  // Get the number of discrete variables
  inline int GetDSplitDim() {
    return dsplitDim;
  }

  // Get domain sizes of discrete variables
  inline const arma::ivec& GetDDomainSize() {
    return dDomainSize;
  }

  /** Decides what branch to to follow.
  @param Dvars    discrete inputs
  @param Cvars    continuous inputs
  @return 0 for left branch, 1 for right branch
  */
  int ChooseBranch(const arma::ivec Dvars, const arma::vec Cvars) {
    return 0;
  }
  /** Computes probability to take first branch. Need this to accomodate probabilistic splitters
  @param Dvars    discrete inputs
  @param Cvars    continuous inputs
  @return probability to take first branch
  */
  double ProbabilityFirstBranch(const arma::ivec Dvars, const arma::vec Cvars) {
    if (ChooseBranch(Dvars, Cvars) == 0)
      return 1;
    else
      return 0;
  }

  // Initializes the data structures used in split variable selection.
  virtual void InitializeSplitStatistics(void) = 0;

  /** Updates the necessary statistics for split variable selection
  @param Dvars    discrete inputs
  @param Cvars    continuous inputs
  @param p1I       probability for the input to have the first class label
  @param p2I       probability for the input to have the second class label
  @param probability  weighting for this datapoint
  */
  virtual void UpdateSplitStatistics(const arma::ivec Dvars, const arma::ivec Cvars,
                                     double p1I, double p2I, double probability) {
  }

  /** Decides on a split variable and frees the data structures used
  in split selection.
  @param type   type of splitting procedure, splitter specific
  @return 0 if a split variable could be computed, -1 otherwise.
  */
  inline int ComputeSplitVariable(int type) {
    return 0;
  }

  /** Is it worth doing more splits in the future?
  @param branch         branch (0 left, 1 right) to which the question applies
  @param Min_no_datapoints  the minimum number of datapoints in order to split
  @return true if can split, false otherwise
  */
  inline bool MoreSplits(int branch, int Min_no_datapoints) {
    return false;
  }

  /** Cleans up all the unnecessary statistics after the decision has been made */
  void DeleteTemporaryStatistics() {
  }

  ~BinarySplitter(void) {
    // dealocate all the resources
  }

  Json::Value ToJson() {
    return Json::Value();
  }
};

#endif // _CLUS_BINARYSPLITTER_H
