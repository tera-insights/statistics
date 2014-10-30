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

// -*- C++ -*-

#if !defined _CLUS_STATISTICSGATHERERS_H
#define _CLUS_STATISTICSGATHERERS_H

#include "general.h"
#include "splitpointcomputation.h" // for gini and split point computation

#include <vector>
#include <math.h>
#include <armadillo>

#ifdef DEBUG_PRINT
#include <iostream>
#endif

namespace CLUS
{

/** Base class for Binomial Statistics. Implements all the
functionality except the computation of the split point.

Such statistics are useful to characterize discrete attribute
variables -- one is constructed for each attribute.
*/
class BasicBinomialStatistics {
 protected:

  // Counts for each value of the attribute for class 0
  arma::vec countsC0;

  // Counts for each value of the attribute (total)
  arma::vec counts;

  // Is the gini gain and split already computed?
  bool gainComputed;

  /** Internal function to compute the split point. Will be
    overloaded in descendants.
    @param N     number of samples used to construct the statistics
    @param alpha_1   apriory probability to have the first class label
    @return best gini gain
   */
  inline double ComputeSplitPoint(double N, double alpha_1) {
    return 0.0;
  }

public:
  /** Default constructor:
    @param DomainSize    size of the domain of the discrete attribute
  */
  BasicBinomialStatistics(int DomainSize = 0)
      : countsC0(DomainSize),
        counts(DomainSize),
        gainComputed(false){
  }

  /** Copy Constructor */
  BasicBinomialStatistics(const BasicBinomialStatistics& o)
      : countsC0(o.countsC0),
        counts(o.counts),
        gainComputed(o.gainComputed) {
  }

  // Get the total nubmer of examples used to build the statistics
  inline double getCount() {
    return arma::sum(counts);
  }

  /** (Re)set the sizse of the domain of the attribute
    @param DomainSize  the new size
  */
  void ResetDomainSize(int DomainSize) {
    countsC0.zeros(DomainSize);
    counts.zeros(DomainSize);
  }

  /** Update the statistics given a data-point
    @param value    value of the attribute in the data-point
    @param classLabel   the class label of the data-point
    @param probability  weitht for the data-point
  */
  void UpdateStatistics(int value, int classLabel, double probability = 1) {
    assert(probability >= 0 && probability <= 1);

    if (classLabel == 0)
      countsC0[value] += probability;
    counts[value] += probability;

    // std::cout << " value=" << value <<  " counts[value]=" << counts[value]
    //           << " countsC0[value]=" << countsC0[value] << std::endl;

    assert(std::isfinite(counts[value]));
  }

  /** Probabilistic version of UpdateStatistics (no class label)
    @param value    value of the attribute in the data-point
    @param p0       probability to have the first class label
    @param p1       probability to have the second class label
  */
  void UpdateStatisticsP(int value, double p0, double p1) {
    countsC0[value] += p0;
    counts[value] += p0 + p1;
    assert(std::isfinite(counts[value]));
  }

  void Merge(const BasicBinomialStatistics& other) {
    countsC0 += other.countsC0;
    counts += other.counts;
  }
};

/** Deterministic version of BasicBinomialStatistics.
In particular, the split is deterministic, i.e. a set.
*/
class BinomialStatistics : public BasicBinomialStatistics {
 protected:
  // the values in the left set; the split set
  std::vector<bool> Split;

  double ComputeSplitPoint(double N, double alpha_1) {
    // use the appropriate gini maximization function
    return DiscreteGiniGain(countsC0, counts, N, alpha_1, Split);
  }

 public:
  /** Default constructor:
    @param DomainSize    size of the domain of the discrete attribute
  */
  BinomialStatistics(int DomainSize = 0)
      : BasicBinomialStatistics(DomainSize),
        Split(0) {
  }

  /** Copy Constructor */
  BinomialStatistics(const BinomialStatistics& o)
      : BasicBinomialStatistics(o),
        Split(o.Split) {
  }

  /** Computes, if necessary and returns the split set.
    @return split set (vector of values in the domain for which left branch is followed)
  */
  std::vector<bool>& GetSplit(void) {
    if (gainComputed)
      return Split;
    else {
      ComputeGiniGain();
      return Split;
    }
  }

  /** Compute the best split and its corresponding gini gain.
    @return best gini gain
  */
  double ComputeGiniGain(void) {
    // create a temporary vector of doubles out of countsC0 so that
    // DiscreteGiniGain function can be used
    gainComputed = true;
    double N = sum(counts); // keeps the total count
    double NC0 = sum(countsC0);

    // std::cout << "Discrete. ComputeGiniGain " << N << " " << NC0 << "\t";

    double alpha_1 = NC0 / N;
    //std::cerr << alpha_1 << "?" << std::endl;
    if (alpha_1 != alpha_1 || N == 0 ) {
      //std::cerr << "no discrete data"<< std::endl;
      return 0;
    }

    double val = DiscreteGiniGain(countsC0, counts, N, alpha_1, Split);

    // std::cout << val << std::endl;

    return val;
  }

};

/** Probabilistic version of BasicBinomialStatistics. In
particular the split is probabilistic.
*/
class ProbabilisticBinomialStatistics : public BasicBinomialStatistics {
 protected:
  // probability to follow the left branch for each possible value
  arma::vec probSet;

 public:
  /** Default constructor:
    @param DomainSize    size of the domain of the discrete attribute
  */
  ProbabilisticBinomialStatistics(int DomainSize = 0)
      : BasicBinomialStatistics(DomainSize),
        probSet(DomainSize) {
  }

  ProbabilisticBinomialStatistics(const ProbabilisticBinomialStatistics& o)
      : BasicBinomialStatistics(o),
        probSet(o.probSet) {
  }

  /** Computes, if necessary, and returns the probability vector.
    @return probability vector  (vector of probabilities to follow left branch for each possible value)
  */
  arma::vec& GetProbabilitySet(void) {
    if (gainComputed)
      return probSet;
    else {
      ComputeGiniGain();
      return probSet;
    }
  }

  /** Compute the best split and its corresponding gini gain.
    @return best gini gain
  */
  double ComputeGiniGain(void) {
    // create a temporary vector of doubles out of countsC0 so that
    // DiscreteGiniGain function can be used
    gainComputed = true;
    double N = sum(counts);
    double NC0 = sum(countsC0);

    // std::cout << "Discrete. ComputeGiniGain " << N << " " << NC0 << "\t";

    double alpha_1=NC0/N;
    //std::cerr << alpha_1 << "?" << std::endl;
    if (alpha_1 != alpha_1 || N== 0) {
      //std::cerr << "no discrete data"<< std::endl;
      return 0;
    }

    double val = ProbabilisticDiscreteGiniGain(countsC0, counts, N, alpha_1, probSet);

    // std::cout << val << std::endl;

    return val;
  }

};

/** Class to maintain statistics for a single continuous
attribute.  Maintains normal approximations to the
distribution of points with each of the two class labels and
uses it to find the best split point and gini gain.
*/
class NormalStatistics {
 protected:
  // number of points with the first class label
  double countC0;

  // number of points with the second class label
  double countC1;

  // sum of points with the first class label
  double sumC0;

  // sum of points with the second class label
  double sumC1;

  // sum squares of points with the first class label
  double sum2C0;

  // sum squares of points with the second class label
  double sum2C1;

  // the computed split point
  double split;

  // solution returned by QDA
  int whichSol;

  // estimate of variance of the split point
  double splitVariance;

  // is the best split determined?
  bool splitVarComputed;

 public:
  /** Default constructor. Must take an int parameter since we
    want to use it as elemet in vectors.
  */
  NormalStatistics(int dummy = 0) {
    Reset();
  }

  /** Copy Constructor */
  NormalStatistics(const NormalStatistics& o)
      : countC0(o.countC0),
        countC1(o.countC1),
        sumC0(o.sumC0),
        sumC1(o.sumC1),
        sum2C0(o.sum2C0),
        sum2C1(o.sum2C1),
        split(o.split),
        whichSol(o.whichSol),
        splitVariance(o.splitVariance),
        splitVarComputed(o.splitVarComputed) {
  }

  // Initializes the values of the statistics
  void Reset(void) {
    countC0 = 0.0;
    countC1 = 0.0;
    sumC0 = 0.0;
    sumC1 = 0.0;
    sum2C0 = 0.0;
    sum2C1 = 0.0;
    splitVarComputed = false;
  }

  // Prints the statistics (debuging)
  void Print(void) {
    std::cout << "Statistics: " << countC0 << " " << countC1 << " ";
    std::cout << sumC0 << " " << sumC1 << " " << sum2C0 << " " << sum2C1 << std::endl;
  }

  // Computes the variance of the combined data (for both classes)
  inline double getVariance(void) {
    return (sum2C0 + sum2C1) / (countC0 + countC1)
        - pow2((sumC0 + sumC1) / (countC0 + countC1));
  }

  // Computes the variance of the split point
  double getSplitVariance(void) {
    if (!splitVarComputed) {
      ComputeGiniGain();
    }
    //return 1.0e-32; // FIXTHIS
    return splitVariance;// is this variance???
  }

  /** Update the statistics given a data-point
    @param value    value of the attribute in the data-point
    @param classLabel   the class label of the data-point
    @param probability  weitht for the data-point
  */
  void UpdateStatistics(double value, int classLabel, double probability = 1.0) {
    if (classLabel==0) {
      countC0 += probability;
      sumC0 += probability * value;
      sum2C0 += probability * pow2(value);
    } else {
      countC1 += probability;
      sumC1 += probability * value;
      sum2C1 += probability * pow2(value);
    }
    assert(std::isfinite(sumC0));// NaN
  }

  void Merge(const NormalStatistics& other) {
    countC0 += other.countC0;
    sumC0 += other.sumC0;
    sum2C0 += other.sum2C0;

    countC1 += other.countC1;
    sumC1 += other.sumC1;
    sum2C0 += other.sum2C0;
  }

  /** Compute the best split and its corresponding gini gain.
    @return best gini gain
  */
  double ComputeGiniGain(void) {
    if (countC0 + countC1 == 0)
      return 0.0;

    double alpha_1 = (double) countC0 / (countC0 + countC1);
    double alpha_2 = (double) countC1 / (countC0 + countC1);
    double eta1 = sumC0 / countC0;
    double eta2 = sumC1 / countC1;
    double var1 = sum2C0 / countC0 - pow2(eta1);
    double var2 = sum2C1 / countC1 - pow2(eta2);

    // std::cout << "ComputeGiniGain " << alpha_1 << " " << alpha_2 << " " << eta1 << " ";
    // std::cout << eta2 << " " << var1 << " " << var2 << "\t";

    split = UnidimensionalQDA(alpha_1, eta1, var1, alpha_2, eta2, var2, whichSol);
    assert(std::isfinite(split));
    splitVariance = UnidimensionalQDAVariance(countC0, eta1, var1,
                                              countC1, eta2, var2, whichSol);
    splitVarComputed = true;

    // compute gini like in function ComputeSeparatingHyperplane_Anova
    double p11 = alpha_1 * PValueNormalDistribution(eta1, sqrt(var1), split);
    double p1_= p11 + alpha_2 * PValueNormalDistribution(eta2, sqrt(var2), split);
    double gini = BinaryGiniGain(p11, alpha_1, p1_);

    // std::cout << split << " " << gini << std::endl;

    return gini;
  }

  // Gets the split point
  inline double GetSplit(void) {
    return split;
  }

};


/** Class implements an approximation of the data with two
  multidimentional normal distributions (Gaussians), one for each
  class label. Useful for finding split points and separating
  hyperplanes.
*/
class MultidimNormalStatistics {
 protected:
  // Dimensionality of the datapoints
  int dim;

  // Sum p1i, Sum p2i
  double s_p1, s_p2;

  // Sum p1i*x, Sum p2i*x
  arma::vec s_p1_x, s_p2_x;

  // Sum p1i*xx^T, Sum p2i*xx^T
  arma::mat S_p1_xxT, S_p2_xxT;

  // The split variable. If -1, oblique split
  int SplitVariable;

  // Equation of the separating hyperplane. Also encodes unidimensional splits
  arma::vec SeparatingHyperplane;

 public:
  // Default constructor
  MultidimNormalStatistics()
      : dim(0),
        s_p1(0),
        s_p2(0),
        s_p1_x(),
        s_p2_x(),
        S_p1_xxT(),
        S_p2_xxT(),
        SplitVariable(0),
        SeparatingHyperplane() {
  } // all the default constructors called for vectors and matrices

  /** Copy Constructor */
  MultidimNormalStatistics(const MultidimNormalStatistics& o)
      : dim(o.dim),
        s_p1(o.s_p1),
        s_p2(o.s_p2),
        s_p1_x(o.s_p1_x),
        s_p2_x(o.s_p2_x),
        S_p1_xxT(o.S_p1_xxT),
        S_p2_xxT(o.S_p2_xxT),
        SplitVariable(o.SplitVariable),
        SeparatingHyperplane(o.SeparatingHyperplane) {
  }

  // Method of computatio of the split point of hyperplane
  enum SeparationType {ANOVA = 0, LDA = 1, QDA = 2};

  // Get the sum of probabilities to belong to the first Gaussian
  double GetS_P1(void) {
    return s_p1;
  }

  // Get the sum of probabilities to belong to the second Gaussian
  double GetS_P2(void) {
    return s_p2;
  }

  /** Changes the dimention of the data and reinitializes everything
    @param Dim  new dimension
  */
  void Resize(int Dim) {
    std::cout << "Size set for cont split stats. Size: " << Dim << std::endl;
    if (dim != Dim) {
      dim = Dim;

      s_p1_x.set_size(dim);
      s_p2_x.set_size(dim);
      S_p1_xxT.set_size(dim, dim);
      S_p2_xxT.set_size(dim, dim);
      SeparatingHyperplane.set_size(dim + 1);
    }

    if (dim > 0) /* otherwise somebody is just freeing the memory */
      Reset();
  }

  // Reset(initialize) the maintained statistics
  void Reset(void) {
    s_p1.fill(0);
    s_p2.fill(0);
    s_p1_x.fill(0);
    s_p2_x.fill(0);
    S_p1_xxT.fill(0);
    S_p2_xxT.fill(0);
  }

  /** Update the statistics maintained with the current datapoint.
    @param values   vector containing the current datapoint
    @param p1     probability for the datapoint to belong to first Gaussian
    @param p2     probability for the datapoint to belong to second Gaussian
  */
  void UpdateStatistics(const arma::vec values, double p1, double p2) {
    if (!finite(p1 + p2))
      return;

    s_p1 += p1;
    s_p2 += p2;

    s_p1_x += p1 * values;
    s_p2_x += p2 * values;

    S_p1_xxT += p1 * values * values.t();
    S_p2_xxT += p2 * values * values.t();
  }

  void Merge(const MultidimNormalStatistics& other) {
    s_p1 += other.s_p1;
    s_p2 += other.s_p2;

    s_p1_x += other.s_p1_x;
    s_p2_x += other.s_p2_x;

    S_p1_xxT += other.S_p1_xxT;
    S_p1_xxT += other.S_p2_xxT;
  }

  /** Changes the statistics to their final form. Must be called
    before the split point can be computed */
  void UpdateParameters() {
    // compute mu1, mu2 in s_p1_x and s_p2_x
    s_p1_x /= s_p1;
    s_p2_x /= s_p2;

    // compute Sigma1 and Sigma2 in S_p1_xxT and S_p2_xxT
    S_p1_xxT = S_p1_xxT / s_p1 - s_p1_x * s_p1_x.t();
    S_p2_xxT = S_p2_xxT / s_p2 - s_p2_x * s_p2_x.t();

    // fix in S_xxT the almost 0.0 entries. This helpes a lot
    S_p1_xxT.transform([&](double val) {
      return (fabs(val) < SMALL_POZ_VALUE) ? 0 : val;
    });
    S_p2_xxT.transform([&](double val) {
      return (fabs(val) < SMALL_POZ_VALUE) ? 0 : val;
    });
  }

  // Get the split variable. If oblique split -1
  int GetSplitVariable(void) {
    return SplitVariable;
  }

  // Get the separating hyperplane (encodes also the split point)
  arma::vec& GetSeparatingHyperplane(void) {
    return SeparatingHyperplane;
  }


  /** Compute the best gini gain for the distribution.
    In the process, the best split point or separating hyperplane is also computed.
    Three types of splits are defined
    LDA:  linear discriminant analysis (fisher discriminant+unidimensional QDA)
    QDA:  linear approximation to quadratic discriminant analysis solution
    ANOVA: best unidimensional split
    @param type  type of split
    @return gini of the split
  */
  double ComputeGiniGain(int type = LDA) {
    double gini, alpha_1, alpha_2;

    if (dim <=1 )
      return 0;

    SplitVariable = -1;

    if (s_p1 + s_p2 == 0)
      return 0.0;

    // compute the best oblique split
    alpha_1 = s_p1 / (s_p1 + s_p2);
    alpha_2 = s_p2 / (s_p1 + s_p2);
    double mass = s_p1 + s_p2;

    UpdateParameters();

    /* All the SeparatingHyperplane computation funcions are normalized
       with respect to variance of the split point */

    switch (type) {
      case ANOVA:
        gini = ComputeSeparatingHyperplane_Anova(
            mass, alpha_1, s_p1_x, S_p1_xxT,
            alpha_2, s_p2_x, S_p2_xxT, SeparatingHyperplane);
        break;
      case LDA:
        gini = ComputeSeparatingHyperplane_LDA(
            mass, alpha_1, s_p1_x, S_p1_xxT,
            alpha_2, s_p2_x, S_p2_xxT, SeparatingHyperplane);
        break;
      case QDA:
        gini = ComputeSeparatingHyperplane_QDA(
            mass, alpha_1, s_p1_x, S_p1_xxT,
            alpha_2, s_p2_x, S_p2_xxT, SeparatingHyperplane);
        break;
      default:
        gini = 0.0; // should never be reached
    }

    SplitVariable = -1;

    arma::uvec indices = find(SeparatingHyperplane.subvec(1, dim) != 0);
    int posSplitVar = (indices.n_elem > 0) ? indices(0) : -1;
    bool isSimple = (indices.n_elem < 2);

    if (isSimple && posSplitVar != -1) {
      SplitVariable = -(posSplitVar + 2);
      std::cout << "Split variable is simple=" << posSplitVar << std::endl;
    } else {
      std::cout << "Simple separation detection failed at: " << posSplitVar << std::endl;
    }

    return gini;
  }

  /** Compute the maximum possible value of the gini
    @return maximum gini
  */
  inline double MaxGini() {
    return 2 * s_p1 * s_p2 / pow2(s_p1 + s_p2);
  }

};

}

#endif // _CLUS_STATISTICSGATHERERS_H
