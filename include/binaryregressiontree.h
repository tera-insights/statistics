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

#if !defined _CLUS_BINARY_REGRESSIONTREE_H_
#define _CLUS_BINARY_REGRESSIONTREE_H_

#include <iostream>
#include <memory>
#include <armadillo>
#include <stdexcept>
#include <cstddef>
#include <jsoncpp/json/json.h>

#include "regressiontree.h"
#include "binaryregressiontreenode.h"

namespace CLUS {

/** Implements a regression tree. The class has three template parameters:
@param T_Distribution     distribution of the elements of the mixture in EM algorithm
@param T_Regressor        regressor model to be used in the leaves. For example constant, linear, etc
@param T_Splitter         takes care of finding the split variable and point once the regression problem is reduced
                          to a classification problem

BinaryRegressionTree takes care of the issues related to the
overall tree, in particular it interfaces with the rest of
the system and preprocesses the data. The actual tree
construction and inference takes place in the nodes of type
BinaryRegressionTreeNode.

*/

template<class T_Distribution, class T_Splitter>
class BinaryRegressionTree : public RegressionTree {
 private:
  using TreeType = BinaryRegressionTree<T_Distribution, T_Splitter>;

  using NodeType = BinaryRegressionTreeNode<T_Distribution, T_Splitter>;
  using NodePtr = std::unique_ptr<NodeType>;

  /// Root of the tree. Specifies the whole tree
  NodePtr root;

  /// list of discrete domain sizes
  DiscreteVector dDomainSize;

  /// number of discrete variables
  SizeType dsplitDim;

  /// number of split variables
  SizeType csplitDim;

  /// number of regression variables
  SizeType regDim;

  /// the minimum number of datapoints in a node to split further
  SizeType min_no_datapoints;

  /// type of split to be passed to splitter, splitter dependent
  int splitType;

  /// threshold for considering a branch irrelevant
  double threshold;

  /// Convergence limit
  double conv_limit;

 public:
  /** Constructor:
    @param DDomainSize   vector containing sizes of domains for discrete variables
    @param CsplitDim     number of split continuous variables
    @param RegDim      number of regressor variables
  */

  BinaryRegressionTree(
      const DiscreteVector& DDomainSize,
      const SizeType CsplitDim,
      const SizeType RegDim,
      const int _splitType = 0,
      const double _threshold = 0.01,
      const double _conv_limit = 1.0e-6,
      const SizeType emRestarts = 3,
      const SizeType emDiscIter = 2,
      const SizeType emLearnIter = 30,
      const SizeType emMaxTrials = 3)
      : root(nullptr),
        dDomainSize(DDomainSize),
        dsplitDim(DDomainSize.n_elem),
        csplitDim(CsplitDim),
        regDim(RegDim),
        min_no_datapoints(10),
        splitType(_splitType),
        threshold(_threshold),
        conv_limit(_conv_limit) {
    using EMConfig = typename NodeType::ConfigType;
    const EMConfig emConfig = {emRestarts, emDiscIter, emLearnIter, emMaxTrials};
    root.reset(new NodeType(1, DDomainSize, CsplitDim, RegDim, emConfig));
  }

  // Copy constructor
  BinaryRegressionTree(const TreeType& other)
      : root(new NodeType(*(other.root))),
        dDomainSize(other.dDomainSize),
        dsplitDim(other.dsplitDim),
        csplitDim(other.csplitDim),
        regDim(other.regDim),
        min_no_datapoints(other.min_no_datapoints),
        splitType(other.splitType),
        threshold(other.threshold),
        conv_limit(other.conv_limit) {
  }

  virtual ~BinaryRegressionTree(void) {}

  // Copy Assignment
  BinaryRegressionTree& operator =(const TreeType& other) {
    root.reset(new NodeType(*(other.root)));
    dDomainSize = other.dDomainSize;
    dsplitDim = other.dsplitDim;
    csplitDim = other.csplitDim;
    regDim = other.regDim;
    min_no_datapoints = other.min_no_datapoints;
    splitType = other.splitType;
    threshold = other.threshold;
    conv_limit = other.conv_limit;

    return *this;
  }

  virtual void StartLearningRound(void) override {
    root->StartLearningEpoch();
  }

  virtual void LearnSample(
      const DiscreteVector& dVars,
      const ContinuousVector& cVars,
      double probability = 1.0) override {
    root->LearnSample(dVars, cVars, probability, threshold);
  }

  virtual void Merge(const TreeType& other) {
    root->Merge(*(other.root));
  }

  virtual bool StopLearningRound(void) override {
    return root->StopLearningEpoch(splitType, min_no_datapoints, conv_limit);
  }

  virtual void StartPruningRound(void) override {
    root->InitializePruningStatistics();
  }

  virtual void PruningSample(
      const DiscreteVector& dVars,
      const ContinuousVector& cVars,
      double y,
      double probability = 1.0) override {
    root->UpdatePruningStatistics(dVars, cVars, y, probability, threshold);
  }

  virtual bool StopPruningRound(void) override {
    return root->FinalizePruningStatistics();
    root->PruneSubtree();
    return true;
  }

  virtual ContinuousValue Infer(
      const DiscreteVector& dVars,
      const ContinuousVector& cVars,
      const double infer_threshold,
      const SizeType maxNodeID = std::numeric_limits<SizeType>::max()) override {
    if (!root) {
      throw std::logic_error("Attempting to infer from empty tree");
    }

    return root->Infer(dVars, cVars, maxNodeID, infer_threshold);
  }

  virtual Json::Value ToJson() const {
    return root->ToJson();
  }

  virtual void ComputeSizesTree(int& nodes, int& term_nodes) const {
    root->ComputeSizesTree(nodes, term_nodes);
  }
};

}

#endif // _CLUS_BINARY_REGRESSIONTREE_H_
