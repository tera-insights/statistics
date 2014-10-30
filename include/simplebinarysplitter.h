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

#if !defined _CLUS_SIMPLEBINARYSPLITTER_H_
#define _CLUS_SIMPLEBINARYSPLITTER_H_

#include "statisticsgatherers.h"
#include <armadillo>

namespace CLUS
{

/** Splitter for traditional classification trees.

    Keeps track of the statistics for variables, decides what the
split predicate should be and indicates split when input is
provided.
*/
class SimpleBinarySplitter
{
protected:

    /// number of discrete variables
    int dsplitDim;

    /// number of continuous variables
    int csplitDim;

    /// vector of discrete domain sizes
    const arma::ivec& dDomainSize;

    /** Indicates on what variable this node splits on. */
    int SplitVariable;

    /** split point if continuous variable is split variable */
    double splitPoint;

    /** The list of values for the left child for the separating variable if
        the split is not oblique. Contains values in order so that
        binary search can be used. */
    std::vector<bool> SeparatingSet;

    /// statistics for discrete variables
    std::vector<BinomialStatistics> discreteStatistics;

    /// statistics for continuous variables
    std::vector<NormalStatistics> continuousStatistics;

    /// number of datapoints
    int count;

    /// number of datapoints with first class label
    int countC0;
public:
    /** Construct object when dimensionality is known
        @param DDomainSize     vector of domain sizes for discrete variables
        @param CsplitDim       number of continuous variables
    */
    SimpleBinarySplitter(const arma::ivec& DDomainSize,int CsplitDim)
        : dsplitDim(DDomainSize.n_elem - 1),
          csplitDim(CsplitDim),
          dDomainSize(DDomainSize) {
    }

    ~SimpleBinarySplitter(void)
    {}

    /// Gets the number of continuous variables
    int GetCSplitDim(void)
    {
        return csplitDim;
    }

    /// Gets the number of discrete variables
    int GetDSplitDim(void)
    {
        return dsplitDim;
    }

    /// Get domain sizes of discrete variables
    const arma::ivec& GetDDomainSize(void)
    {
        return dDomainSize;
    }

    /** Initializes the datastructures used in split variable selection */
    void InitializeSplitStatistics(void)
    {
        count = 0;
        countC0 = 0;

        discreteStatistics.resize(dsplitDim);
        for (int i = 0; i < dsplitDim; i++)
        {
            discreteStatistics[i].ResetDomainSize(dDomainSize[i]);
        }

        continuousStatistics.resize(csplitDim);
    }


    /** Decides what branch to to follow.
        @param Dvars        discrete inputs
        @param Cvars        continuous inputs
        @return 0 for left branch, 1 for right branch
    */
    int ChooseBranch( const arma::ivec& Dvars, const arma::vec& Cvars)
    {
        if (SplitVariable <=-1)
        {
            // split on a continuous variable
            // the split variable is actually -SplitVariable-1
            if (Cvars[-SplitVariable-1]<=splitPoint)
                return 0;
            else
                return 1;
        }
        else
        {
            // split on a discrete variable
            int value=Dvars[SplitVariable];
            // look for value in SeparatingSet
            if (IsPointInSet<bool>((bool) value, SeparatingSet))
                return 0;
            else
                return 1;
        }
    }

    /** Updates the necessary statistics for split variable selection
        @param Dvars        discrete inputs
        @param Cvars        continuous inputs
        @param classLabel   class label of the datapoint
    */
    void UpdateSplitStatistics( const arma::ivec& Dvars, const arma::vec& Cvars,
                                int classLabel)
    {
#ifdef DEBUG_PRINT
        cout << "CL=" << classLabel << " ";
#endif

        count++;
        if (classLabel == 0)
            countC0++;

#ifdef DEBUG_PRINT

        std::cout << "count=" << count << " countC0=" << countC0 << std::endl;
#endif

        for (int i=0; i<dsplitDim; i++)
        {
            int value=Dvars[i];
            discreteStatistics[i].UpdateStatistics(value, classLabel);
        }

        // update continuous statistics
        for (int i=0; i<csplitDim; i++)
        {
            double value=Cvars[i];
            continuousStatistics[i].UpdateStatistics(value, classLabel);
        }
    }

    /** Decides on a split variable and frees the data structures used
        in split selection.
        @return true if a split variable could be computed, false otherwise.
    */
    bool ComputeSplitVariable(void)
    {
        // make node a leaf if not enough data to take a split decision
        if (countC0<2.0 || count-countC0<2.0)
        {
#ifdef DEBUG_PRINT
            cout << "Making the node a leaf. countC0=" << countC0
            << " count=" << count << endl;
#endif

            return false;
        }

        double maxgini=0.0;

        // go over the discrete attributes and find the best one
        for (int i=0; i<dsplitDim; i++)
        {
            double curr_gini=discreteStatistics[i].ComputeGiniGain();
#ifdef DEBUG_PRINT

            cout << "\tVariable: " << i << " gini=" << curr_gini << endl;
#endif

            if (curr_gini>maxgini)
            {
                maxgini=curr_gini;
                SplitVariable=i;
            }
        }

        // go over continuous variable
        for (int i=0; i<csplitDim; i++)
        {
            double curr_gini=continuousStatistics[i].ComputeGiniGain();
#ifdef DEBUG_PRINT

            cout << "\tVariable: " << (-i-1) << " gini=" << curr_gini << endl;
#endif

            if (curr_gini>maxgini)
            {
                maxgini=curr_gini;
                SplitVariable=-(i+1);
            }
        }

        if (maxgini==0.0)
            return false; // make the node a leaf, nobody can do a reasonable split

        // not set the split point for the variable picked as the split point
        if (SplitVariable >= 0)
            SeparatingSet=discreteStatistics[SplitVariable].GetSplit();
        else
        {
            int splitVar=-SplitVariable-1;
            splitPoint=continuousStatistics[splitVar].GetSplit();
        }

#ifdef DEBUG_PRINT
        cout << "Split variable is " << SplitVariable << " and split point ";
        if (SplitVariable>=0)
            cout << SeparatingSet << endl;
        else
            cout << splitPoint << endl;
#endif

        return true;
    }

    /// Computes the majority class label, 0 or 1
    int ComputeClassLabel(void)
    {
        if (countC0>=count-countC0)
            return 0;
        else
            return 1;
    }

    /** Is it worth doing more splits in the future?
        @param minMass    the minimum number of datapoints in order to split
        @param nodeId     the ID of the node that owns the splitter
        @return true if can split, false otherwise
    */
    bool MoreSplits(int minMass, int nodeId)
    {
        return ( count>=minMass && countC0!=0 && countC0!=count );
    }

    /** Save the state of the splitter to an output stream
        @param out    output stream
    */
    void SaveToStream(std::ostream& out)
    {
        out << SplitVariable << " ";
        if (SplitVariable>=0)
        {
            out << "( " << count << " ) [ ";
            for (int i=0; i<SeparatingSet.size(); i++)
                out << SeparatingSet[i] << " ";
            out << " ] ";
        }
        else
        {
            out << "( " << splitPoint << " ) ";
        }
    }
};
}

#endif // _CLUS_SIMPLEBINARYSPLITTER_H_
