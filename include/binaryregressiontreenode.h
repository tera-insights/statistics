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

#ifndef _CLUS_REGRESSIONTREENODE_H_
#define _CLUS_REGRESSIONTREENODE_H_

#include <memory>
#include <armadillo>
#include <array>
#include <cstddef>
#include <cstdint>
#include <jsoncpp/json/json.h>

#include "distribution.h"
#include "em.h"
#include "regressor.h"

namespace CLUS {
/** Class used in building regression trees. Supports both classic and
    probabilistic regression trees through different splitters.
    Every class is a node in the tree, there can be 0 or 2 children
    NodeId determines the position of the node in the tree. The nodeId
    od children is alwais 2*nodeId and 2*nodeId+1. The nodeId starts from 1
    for the root.

    @param T_Distribution     the distribution that approximates the data
    @param T_Regressor        the regressor that commes with T_Distribution
    @param T_Spliter          the class that in a final scan can compute the splitting predicate
*/
template<class T_Distribution, class T_Splitter>
class BinaryRegressionTreeNode {
 public:
  using SizeType = std::size_t;

  using ExpectMaxer = CLUS::ExpectationMaximizer<T_Distribution>;

  using ContinuousValue = double;
  using ContinuousVector = arma::vec;

  using DiscreteValue = arma::sword;
  using DiscreteVector = arma::ivec;

  using ConfigType = typename ExpectMaxer::ConfigType;

 private:
  using NodeType = BinaryRegressionTreeNode<T_Distribution, T_Splitter>;
  using NodePtr = std::unique_ptr<NodeType>;
  using ChildArray = std::array<NodePtr, 2>;

  using RegressionPtr = std::unique_ptr<Regressor>;
  using ExpectMaxerPtr = std::unique_ptr<ExpectMaxer>;

 private:
  enum class State {
    em,
    split,
    regression,
    stable
  };

  // unique identifier of the cluster for a regression tree.
  // the id of the children is always 2*Id, 2*Id+1
  SizeType nodeId;

  // how many times EM was attempted
  SizeType emTrials;

  // the number of continuous split variables
  SizeType csDim;

  // the number of regressors
  SizeType regDim;

  // the state of the node. At creation em. At load stable
  State state;

  // the children of this node
  ChildArray Children;

  // splitter for split criterion
  T_Splitter splitter;

  // Overall distribution of the data, used for building the regressor
  T_Distribution overallDist;

  // regressor for the node
  RegressionPtr regressor;

  // Expectation Maximization engine
  ConfigType expectMaxerConfig;
  ExpectMaxerPtr expectMaxer;

  // sum of squared differences between prediction and true value
  double pruningCost;

  // number of samples in this node for pruning
  int pruningSamples;

 public:
  // Constructor for a leaf
  //
  // @param NodeId        ID of this node
  // @param CsDim         number of split variables
  // @param RegDim        number of regressor variables
  BinaryRegressionTreeNode(SizeType NodeId, SizeType CsDim, SizeType RegDim)
      : nodeId(NodeId),
        emTrials(0),
        csDim(CsDim),
        regDim(RegDim),
        state(State::regression),
        Children(),
        splitter(),
        overallDist(RegDim),
        regressor(nullptr),
        expectMaxerConfig(),
        expectMaxer(nullptr),
        pruningCost(0),
        pruningSamples(0) {
  }

  // Constructor for an intermediate node
  //
  // @param NodeId        ID of this node
  // @param DDomainSize     vector of domain sizes of discrete variables
  // @param CsplitDim       number of split variables
  // @param RegDim        number of regressor variables
  BinaryRegressionTreeNode(SizeType NodeId, const DiscreteVector& DDomainSize,
                           SizeType CsplitDim, SizeType RegDim,
                           const ConfigType emConfig)
      : nodeId(NodeId),
        emTrials(0),
        csDim(CsplitDim),
        regDim(RegDim),
        state(State::em),
        Children(),
        splitter(DDomainSize, CsplitDim, RegDim),
        overallDist(RegDim),
        regressor(nullptr),
        expectMaxerConfig(emConfig),
        expectMaxer(new ExpectMaxer(RegDim, emConfig)),
        pruningCost(0),
        pruningSamples(0) {
  }

  // Copy Constructor
  BinaryRegressionTreeNode(const NodeType& other)
      : nodeId(other.nodeId),
        emTrials(other.emTrials),
        csDim(other.csDim),
        regDim(other.regDim),
        state(other.state),
        Children(),
        splitter(other.splitter),
        overallDist(other.overallDist),
        regressor(),
        expectMaxerConfig(other.expectMaxerConfig),
        expectMaxer(),
        pruningCost(other.pruningCost),
        pruningSamples(other.pruningSamples) {
    if (!other.IsLeaf()) {
      Children[0].reset(new NodeType(*(other.Children[0])));
      Children[1].reset(new NodeType(*(other.Children[1])));
    }

    if (other.regressor)
      regressor.reset(other.regressor->clone());

    if (other.expectMaxer)
      expectMaxer.reset(new ExpectMaxer(*(other.expectMaxer)));
  }

  // Destroys recursively the tree
  ~BinaryRegressionTreeNode() {
  }

  // Get the ID of this node
  SizeType GetNodeId() const {
    return nodeId;
  }

 private:
  bool IsLeaf() const {
    for (const auto& child : Children)
      if (child)
        return false;
    return true;
  }

 public:
  // Compute the size of the tree
  //
  // @param nodes     output variable where total number of nodes is placed
  // @param term_nodes  output variable where total number of terminal nodes is placed
  void ComputeSizesTree(int& nodes, int& term_nodes) const {
    nodes++;
    if (!IsLeaf()) {
      Children[0]->ComputeSizesTree(nodes, term_nodes);
      Children[1]->ComputeSizesTree(nodes, term_nodes);
    } else {
      term_nodes++;
    }
  }

  // Begin the learning process. Initializes splitter if necessary
  void StartLearningEpoch() {
    // std::cout << "Starting state " << static_cast<int>(state) << std::endl;
    switch (state) {
      case State::stable:
        if (!IsLeaf()) {
          Children[0]->StartLearningEpoch();
          Children[1]->StartLearningEpoch();
        }
        break;
      case State::em:
        expectMaxer->start_round();
        break;
      case State::split:
        splitter.InitializeSplitStatistics();
        break;
      case State::regression:
        // nothing to do
        break;
    }
  }

  // Learn a data sample, i.e. pass it down the tree or to the splitter
  //
  // @param Dvars     input discrete variables
  // @param Cvars     input continuous variables
  // @param probability   weight of the datapoint
  // @param threshold   minimum probability to follow a branch
  void LearnSample(const DiscreteVector& Dvars, const ContinuousVector& Cvars,
                   double probability, double threshold = .01) {
    double p0, p1, probLeft, probRight;

    if (probability < threshold)
      return;

    // Regression-only continuous variables
    auto regVars = Cvars.subvec(csDim, Cvars.n_elem - 1);

    switch (state) {
      case State::stable:
        // Pass the sample to the right child
        if (IsLeaf())
          return;

        // propagate the learning
        probLeft = probability * splitter.ProbabilityLeft(Dvars, Cvars);
        probRight = probability - probLeft;

        if (probLeft >= threshold)
          Children[0]->LearnSample(Dvars, Cvars, probLeft, threshold);

        if (probRight >= threshold)
          Children[1]->LearnSample(Dvars, Cvars, probRight, threshold);

        break;
      case State::em:
        expectMaxer->learn(regVars, probability);
        break;
      case State::split:
        p0 = expectMaxer->prob_left(regVars, probability);
        p1 = expectMaxer->prob_right(regVars, probability);

        if (p0 + p1 > 0)
          splitter.UpdateSplitStatistics(Dvars, Cvars, p0 / (p0 + p1),
                                         p1 / (p0 + p1), probability);
        break;
      case State::regression:
        overallDist.Update(regVars, probability);
        break;
    }
  }

  // Used for merging two binary regression trees node-by-node.
  void Merge(const NodeType& other) {
    switch (state) {
      case State::stable:
        if (IsLeaf()) {
          Children[0]->Merge(*(other.Children[0]));
          Children[1]->Merge(*(other.Children[1]));
          break;
        }
      case State::em:
        expectMaxer->Merge(*(other.expectMaxer));
        break;
      case State::split:
        splitter.Merge(other.splitter);
        break;
      case State::regression:
        overallDist.Merge(other.overallDist);
        break;
    }
  }

  // Stop the learning process
  //
  // @return        true if stopped successfully
  //            false if still processing data
  bool StopLearningEpoch(int splitType, int min_no_datapoints, double convLimit) {
    bool shouldStop = true;

    switch (state) {
      case State::em:
        expectMaxer->end_round(convLimit);

        if (expectMaxer->done()) {
          if (expectMaxer->should_split()) {
            std::cout << "Decided to split." << std::endl;
            // Have 2 candidate child distributions, try to split on them
            state = State::split;
          } else {
            std::cout << "Decided not to split." << std::endl;
            // No valid distributions, just become a leaf
            state = State::regression;
          }
        }

        shouldStop = false;
        break;
      case State::split:
        // cout << "Node: " << nodeId << " being split" << endl;

        state = State::regression;
        shouldStop = false;

        if (splitter.ComputeSplitVariable(splitType)) {
          for (SizeType i = 0; i < Children.size(); i++) {
            SizeType newID = (nodeId * Children.size()) + i;
            if (splitter.MoreSplits(i, min_no_datapoints)) {
              // Create possible intermediate node
              Children[i].reset(new NodeType(
                  newID,
                  splitter.GetDDomainSize(),
                  splitter.GetCSplitDim(),
                  splitter.GetRegDim(),
                  expectMaxerConfig
              ));
            } else {
              // Create leaf node
              Children[i].reset(new NodeType(
                  newID,
                  splitter.GetCSplitDim(),
                  splitter.GetRegDim()
              ));
            }
          }
        }

        splitter.DeleteTemporaryStatistics();

        shouldStop = false;
        break;
      case State::regression:
        overallDist.Estimate();
        regressor.reset(overallDist.CreateRegressor());
        std::cout << "Created regressor" << std::endl;
        std::cout << "Regressor pointer: " << regressor.get() << std::endl;
        std::cout << "Regressor: " << std::endl << regressor->ToJson() << std::endl;

        state = State::stable;
        shouldStop = IsLeaf();

        break;
      case State::stable:
        if (!IsLeaf()) {
          shouldStop = shouldStop & Children[0]->StopLearningEpoch(
              splitType, min_no_datapoints, convLimit);

          shouldStop = shouldStop & Children[1]->StopLearningEpoch(
              splitType, min_no_datapoints, convLimit);
        }
        break;
    }

    return shouldStop;
  }

  // Initialize statistics about pruning
  void InitializePruningStatistics() {
    pruningCost = 0.0;
    pruningSamples = 0;
    if (!IsLeaf()) {
      Children[0]->InitializePruningStatistics();
      Children[1]->InitializePruningStatistics();
    }
  }

  // Update pruning stats with new data
  //
  // @param Dvars     discrete variables
  // @param Cvars     continuous variables
  // @param y       true output
  // @param probability   weight of the datapoint
  // @param threshold   minimum probability to follow a branch
  void UpdatePruningStatistics(const DiscreteVector& Dvars,
                               const ContinuousVector& Cvars, double y,
                               double probability, double threshold) {
    // pruningTotalMass += probability;

    double predY = regressor->Predict(Cvars + csDim);
    pruningCost += (y - predY) * (y - predY) * probability;

    // update pruning statistics for the proper children
    if (IsLeaf())
      return; // stop the process

    double probLeft = probability * splitter.ProbabilityLeft(Dvars, Cvars);
    // pruningTotalMassLeft += probabilityLeft;
    double probRight = probability - probLeft;

    if (probLeft >= threshold)
      Children[0]->UpdatePruningStatistics(Dvars, Cvars, y, probLeft,
                                           threshold);

    if (probRight >= threshold)
      Children[1]->UpdatePruningStatistics(Dvars, Cvars, y, probRight,
                                           threshold);
  }

  // Postprocessing for prunning statistics
  bool FinalizePruningStatistics() {
    return false;
  }

  // Return the optimal cost for this subtree and cuts the subtree to optimal size
  // @return optimal cost
  double PruneSubtree() {
    // double alpha /* alpha is the cost for a leaf */){
    if (IsLeaf()) {
      // node is a leaf
      return pruningCost;
    } else {
      // node is an intermediary node
      double pruningCostChildren = Children[0]->PruneSubtree()
                                 + Children[1]->PruneSubtree();

      if (pruningCost <= pruningCostChildren) {
        // prune the tree here
        Children[0].reset();
        Children[1].reset();

        return pruningCost;
      } else {
        // tree is good as it is
        return pruningCostChildren;
      }
    }
  }

  // Does the inference
  //
  // @param Dvars     input discrete variables
  // @param Cvars     input continuous variables
  // @param maxNodeId   maximum ID of a node to propagate recursive call (limits size)
  // @param threshold   minimum probability to follow a branch
  double Infer(const DiscreteVector& Dvars, const ContinuousVector& Cvars,
               SizeType maxNodeId, double threshold) {
    if (IsLeaf() || nodeId > maxNodeId) {
      // leaf node or level cut
      return regressor->Predict(Cvars + csDim);
      //return nodeId;
    } else {
      double pChild1 = splitter.ProbabilityLeft(Dvars, Cvars);

      return pChild1 * ((pChild1 >= threshold)
                 ? Children[0]->Infer(Dvars, Cvars, maxNodeId, threshold)
                 : 0)
           + (1 - pChild1) * ((1 - pChild1 >= threshold)
                 ? Children[1]->Infer(Dvars, Cvars, maxNodeId, threshold)
                 : 0);
    }
  }

  Json::Value ToJson() const {
    Json::Value ret;
    if (!IsLeaf()) {
      ret["splitter"] = splitter.ToJson();
      ret["left"] = Children[0]->ToJson();
      ret["right"] = Children[1]->ToJson();
    }
    ret["regressor"] = regressor->ToJson();
    return ret;
  }
};
}

#endif //  _CLUS_REGRESSIONTREENODE_H_
