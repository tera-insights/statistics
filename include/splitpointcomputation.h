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

/* Implement functions for split point computation together with
   necessary auxiliary functions
*/

#if !defined _CLUS_SPLITPOINTCOMPUTATION_H_
#define _CLUS_SPLITPOINTCOMPUTATION_H_

#include <math.h>
#include <iostream>
#include <armadillo>
#include <cassert>
#include <stdlib.h>
#include <boost/math/special_functions/beta.hpp>

#include "general.h"
#include "math.h"

namespace CLUS {

/** Determines in log time if a point is in a set.
@param value   the point
    @param set     the vector of sorted values
@return true if point value in set, false otherwise
*/
template<class T>
bool IsPointInSet(T value, arma::Col<T> set)
{
    int l = 0, r = set.n_elem - 1;

    assert(r >= 0);

    while (l <= r)
    {
        int m = (l + r) / 2;
        T vm = set[m];

        if (vm == value)
            return true;

        if (vm < value)
            l = m + 1;
        else
            r = m - 1;
    }
    return false;
}

/** Compute P[X>=val] for X~Binomial(N,p).
    This is exactly the normalized incomplete beta function (see Eric's encyclopedia).
    We use gsl to get this function. The result is IB_p(val,N-val-1.0)
    @param N    number of coinflips
    @param p    probability to see head
    @param val  value for which the p-value is computed
    @return the probability that at least val heads are observed out of N coin flips
*/
double PValueBinomialDistribution(double N, double p, double val)
{

    double mean = N * p;
    double std = std::sqrt(N * p * (1-p)); // standard deviantion

    // std::cout << "mean=" << mean << " std=" << std << std::endl;

    // if val is not withing 10 standard deviations do not bother with the approximation
    if (val < mean - 10 * std)
        return 1.0;

    if (val > mean + 10 * std)
        return 0.0;

    double result = boost::math::ibeta(val + 1, N - val, p);
    return result;
}

/** Computes int_{x>=eta} N(mu,sigma^2) dx, i.e. the probability a
    sample form the normal distribution with mean mu and standard
    deviation sigma is as large or larger than eta.
*/
double PValueNormalDistribution(double mu, double sigma, double eta)
{
    if (sigma == 0.0){
        if (mu > eta) {
            return 1.0;
        } else {
            if (mu == eta)
                return .5;
            else
                return 0.0;
        }
    }
    return .5 * (1.0 + erf((mu - eta) / sqrt(2 * pow2(sigma))));
}

/** Computes int_{n'*(x-xc)>=0} N(mu,Sigma) dx, i.e. probability that
    a sample from the multidimensional Gaussian distribution is above
    a hyperplane (n,xc).
    @param mu          the mean of the distribution
    @param cholSigma   the Cholesky decomposition of the covariance matrix of the distribution
    @param n           the normal to the hyperplane
    @param xc          a point contained in the hyperplane
    @return the probability a sample from the distribution is above the hyperplane
 */
double PValueNormalDistribution(const arma::vec mu,
                                const arma::mat cholSigma,
                                const arma::vec n,
                                const arma::vec xc)
{
    /* Crafted after pvalue.m and thing.m (Matlab prototyping)
       The final formula is:
       1                   /       / sigma v \ \
              P=------------------------------------- | 1+ Erf| ------- | |
       2 sigma det(cholSigma) Sqrt[ det(S) ] \       \ Sqrt[2] / /

       where:
              v=n'*(mu-xc)/norm(n),
              sigma=Sqrt[ s-w'S^-1w ],
       / s w'\
       (M'*cholSigma*cholSigma'*M)^-1=|     |
       \ w S /
       and M is picked such that  M'n=e1

    */

    int d = n.n_elem; // dimentionality of the space

    // compute v
    double v = arma::dot(n, mu - xc) / sqrt(arma::dot(n, n));

    // build initial M
    arma::mat M(d , d); // already initialized to 0
    /* pick robustly the other d-1 vectors, eliminate ei, i=max(abs(n))
       since is the most "paralel" with n */
    // find i first
    arma::uword i;
    arma::vec abs_n = arma::abs(n);
    abs_n.max(i);
    M.diag().fill(1.0);
    M.col(0) = n;
    // put e1 in row i
    if (i != 0) {
        M(i, i)=0.0;
        M(0, i)=1.0;
    }

    // Gram-Schmidt ortogonalize M
    for (int j = 0; j < d; j++) {
        for (int k = 0; k < j; k++)
            // Decrease row j by (row k) * (row k . row j)
            M.row(j) -= arma::dot(M.row(j), M.row(k)) * M.row(k);
        // Normalize row j
        M.row(j) = normalise(M.row(j));
    }

    // compute cholSigma^-1 M into M
    M = cholSigma.i() * M;

    // compute s
    double s = sum(square(M.diag()));

    // compute w
    arma::vec w  = arma::trans(M.col(0).t() * M.cols(1, d - 1));

    // compute S
    arma::mat S = arma::chol(M.cols(1, d - 1).t() * M.cols(1, d - 1));

    // compute cholS^-1 * w
    arma::vec aux = S.i() * w;

    double sigma = s - arma::dot(aux, aux);

    sigma = sqrt(sigma);
    double det_cholS = arma::det(S);
    double det_cholSigma = arma::det(cholSigma);

    // we can compute now the PValue
    return 1 / (2 * sigma * det_cholSigma * det_cholS) * (1 + erf(sigma * v / sqrt(2.0)));
}

/** Auxiliary type. Used for sorting elements in a discrete domain */
struct T_array
{
    int index;
    double crit;
};

/** Auxiliary function for sort elements. Provides comparition between T_array values.
    @param x      first value
    @param y      second value
    @return true if x<y, false otherwise
*/
int compare_array_elements(const void* x, const void* y)
{
    const T_array ax = *( (T_array*)x );
    const T_array ay = *( (T_array*)y );
    if (ax.crit<ay.crit)
        return -1;
    else
        if (ax.crit == ay.crit)
            return 0;
        else
            return 1;
}

/** Computes gini gain.
    @param p11     probability to have the first class label and be at the left of split
    @param p_1     probability to have the first class label
    @param p1_     probability to be at the left of split
*/
double BinaryGiniGain(double p11, double p_1, double p1_)
{
    return p1_>0.0 ? 2*pow2(p11-p_1*p1_)/( p1_*(1-p1_) ) : 0.0;
}

/** Computes the maximum gain in gini by splitting on a discrete variable and the actual split
    Split: return the best split here
    Return: the new gini

    If d_s_p1[i]=d_N[i]=0 we don't know anything about value i. To avoid biases we distribute
    this values among the two splits.

    Theorem 9.4 from Breiman et al. justifies the linear algorithm.
    @param d_s_p1     vector with the sum of probabilities to have the first class label for each possible value
    @param d_N        vector with the number of samples for each possible value
    @param N          total number of samples
    @param alpha_1    apriory probability to belong to have the first class label
    @param Split      vector where the best split set is stored (returned)
    @return best gini gain
*/
double DiscreteGiniGain(arma::vec& d_s_p1, arma::vec& d_N, double N, double alpha_1,
                        std::vector<bool>& Split)
{

    int n = d_N.n_elem;

    T_array* B_order=new T_array[n];

    for (int i=0; i<n; i++)
    {
        B_order[i].index=i;
        B_order[i].crit=d_N[i]>0 ? d_s_p1[i]/d_N[i] : -1.0;
        ;
    }

    qsort((void*)B_order, n, sizeof(T_array), compare_array_elements);

    //    std::cout << "B_order ";
    //     for (int i=0; i<n; i++)
    //       std::cout << " [" << B_order[i].index << "," << B_order[i].crit << "] ";
    //     std::cout << std::endl;


    // find the splitting point
    int pos=-1;  // the splitting point is between pos and pos+1
    double p11=0.0;
    double p1_=0.0;
    double maxgini=0.0;
    bool partonefirst=true; // if true partition one is at the left
    for (int i=0; i<n-1; i++)
    {
        p11+=d_s_p1[B_order[i].index]/N;
        p1_ +=1.0*d_N[B_order[i].index]/N;
        double gini=BinaryGiniGain(p11,alpha_1,p1_);

        //std::cout << "XXXXXX " << p11 << " " << p1_ << " " << gini << std::endl;

        if (p1_>=0.0 && p1_<=1.0 && gini>=maxgini)
        {
            maxgini=gini;
            pos=i;
            if ( p11<(p1_-p11) )
                partonefirst=false;
            //std::cout << "maxgini=" << maxgini << " i=" << i << std::endl;

        }
    }

    double boundary;
    if (pos==-1 || B_order[pos].crit==-1.0)
        return 0.0;
    else
        boundary=d_s_p1[B_order[pos].index]/d_N[B_order[pos].index];

    // std::cout << "boundary=" << boundary << " pos=" << pos << std::endl;

    int nr=0;
    int buff[MAX_CATEGORIES];
    for (int i=0; i<n; i++) {
        bool add_val = false;

        double p1_o_n = d_s_p1[i] / d_N[i];

        if( d_N[i] == 0.0 ) {
            // Decide by coinflip
            add_val = (rand() * (1.0 / RAND_MAX)) > 0.5;
        } else if( !partonefirst ) {
            // If not part one first, add it if d_N[i] / d_s_p1[i] > boundary
            add_val = is_greater_eps(p1_o_n, boundary);
        } else {
            // If part one first, add it if d_N[i] / d_s_p1[i] <= boundary
            add_val = is_less_equal_eps(p1_o_n, boundary);
        }

        if ( add_val )
        {
            // std::cout << "Adding " << i << " " << nr << std::endl;
            buff[nr]=i;
            nr++;
        }
    }

    //std::cout << "BUFFER ";
    //for (int i=0; i<nr; i++)
    //std::cout << buff[i] << " ";
    //std::cout << std::endl;

    // put the result in Split
    Split.resize(nr);
    for (int i=0; i<nr; i++)
        Split[i] = (bool) buff[i];

    //std::cout << "Split ";
    // for (int i=0; i<Split.size(); i++)
    //std::cout << Split[i] << " ";
    //std::cout << " ----- " << nr << std::endl;

    return (nr>0 && nr<=n-1) ? maxgini : 0.0; // to avoid pseudo nonzero ginis
}

/** Sister function of DiscreteGiniGain. Instead of finding a split set it finds
    probabilities that a point belongs to the left set.
    @param d_s_p1     vector with the sum of probabilities to have the first class label for each possible value
    @param d_N        vector with the number of samples for each possible value
    @param N          total number of samples
    @param alpha_1    apriory probability to belong to have the first class label
    @param probSet    vector where the best probability set is stored (returned)
    @return best gini gain
*/
double ProbabilisticDiscreteGiniGain(const arma::vec& d_s_p1, const arma::vec& d_N,
                                     double N, double alpha_1, arma::vec& probSet) {
    int n = d_N.n_elem;

    T_array* B_order=new T_array[n];

    for (int i = 0; i < n; i++)
    {
        B_order[i].index = i;
        B_order[i].crit = (d_N[i] > 0.0) ? d_s_p1[i] / d_N[i] : -1.0;
    }

    qsort((void*) B_order, n, sizeof(T_array), compare_array_elements);

    // find the splitting point
    int pos = -1;  // the splitting point is between pos and pos+1
    double p11 = 0.0;
    double p1_ = 0.0;
    double maxgini = 0.0;
    bool partonefirst = true; // if true partition one is at the left
    for (int i = 0; i < n-1; i++) {
        p11 += d_s_p1[B_order[i].index] / N;
        p1_ += 1.0 * d_N[B_order[i].index] / N;
        double gini = BinaryGiniGain(p11, alpha_1, p1_);

        if (p1_ >= 0.0 && p1_ <= 1.0 && gini >= maxgini) {
            maxgini = gini;
            pos = i;
            if (p11 < (p1_ - p11))
                partonefirst = false;
        }
    }

    // take the boundary to be the average of the ratios where maximum occurs
    double boundary;
    probSet.set_size(n);
    if (pos == -1 || B_order[pos].crit == -1.0)
    {
        for (int i = 0; i < n; i++)
            probSet[i] = 0.5;
        return 0.0;
    }
    else
        boundary = (d_s_p1[B_order[pos].index]/d_N[B_order[pos].index] +
                    d_s_p1[B_order[pos+1].index]/d_N[B_order[pos+1].index]) / 2.0;

    delete [] B_order;

    // for every element in the domain compute the probability that the
    // observed ratio is at the left of the boundary point assuming Binomial distribution
    for (int i = 0; i < n; i++) {
        if (d_N[i] == 0.0) {
            probSet[i] = .5;
        } else {
            double p = d_s_p1[i] / d_N[i];
            probSet[i] = 1.0 - PValueBinomialDistribution(d_N[i], p, boundary * d_N[i]);
            //cout << "i:" << i << " s_p1=" << d_s_p1[i] << " d_N=" << d_N[i]
            //     << " probSet[]=" << probSet[i] << endl;

        }
    }

    return maxgini;
}

/** Implements unidimensional Quadratic Discriminant Analysis,
    i.e. finds separator between tow unidimensional normal
    distributrions.

    Forms equation eta^2(1/var1 -
    1/var2)-2eta(eta1/var1-eta2/var2)+eta1^2/var1-eta2^2/var2 =
    2ln(alpha1/alpha2)-ln(var1/var2) and solves it. In degenerate
    cases, a linear version of the equation is solved.

    @param alpha_1    probability to belong to the first distribution
    @param eta1       mean of first distribution
    @param var1       variance of first distribution
    @param alpha_2    probability to belong to the second distribution
    @param eta2       mean of second distribution
    @param var2       variance of second distribution
    @param whichSol   set to 0 if first order equation solved, 1 if first sol of second
                      order equation and 2 for second solution. Is set to 3 if the
		      failover solution is returned
    @return           the best split point
*/
double UnidimensionalQDA(double alpha_1, double eta1, double var1,
                         double alpha_2, double eta2, double var2, int& whichSol)
{
    whichSol = 4; // special solution: eta1 or eta2

    if (std::isnan(alpha_1+eta1+var1))
    {
        assert(!std::isnan(alpha_2 + eta2 + var2));
        return eta2;
    }

    if (std::isnan(alpha_2+eta2+var2))
    {
        assert(!std::isnan(alpha_1 + eta1 + var1));
        return eta1;
    }

    if (!NonEqual(eta1,eta2))
        return eta1;

    double a=1.0/var1-1.0/var2;
    double b=-2.0*(eta1/var1-eta2/var2);
    double c=pow2(eta1)/var1-pow2(eta2)/var2-2.0*log(alpha_1/alpha_2)/log(exp(1.0))+
             log(var1/var2)/log(exp(1.0));


    if (var1<=SMALL_POZ_VALUE || var2<=SMALL_POZ_VALUE)
    {
        whichSol=5;
        return alpha_1*eta1+alpha_2*eta2;
    }

    goto default_solution;

    if ( fabs(a) < 1e-8*(1.0/var1+1.0/var2) )
    {
        // first order equation

        if (NonZero(b))
        {
            double eta_l=-c/b;
            if ( (eta1-eta_l)*(eta2-eta_l)<0.0 )
            {
                whichSol=0;
                return eta_l;
            }
        } // otherwise return default solution; fallover
    }
    else
    {
        // second order equation
        double disc=pow2(b)-4.0*a*c;
        if (disc<0.0) // no solution, return default
            goto default_solution;
        else
        {
            // the two solutions
            double eta_1=(-b+sqrt(disc))/(2.0*a);
            double eta_2=(-b-sqrt(disc))/(2.0*a);
            // find the one in between
            if ( (eta1-eta_1)*(eta2-eta_1) < 0.0)
            {
                whichSol=1;
                return eta_1;
            }
            else
                if ( (eta1-eta_2)*(eta2-eta_2) < 0.0)
                {
                    whichSol=2;
                    return eta_2;
                } // otherwise return default solution; fallover
        }
    }

default_solution:
    whichSol=3; // default solution
    return (alpha_1*sqrt(var1)*eta1+alpha_2*eta2*sqrt(var2))/
           (alpha_1*sqrt(var1)+alpha_2*sqrt(var2));
    //return alpha_1*eta1+alpha_2*eta2;
}

/** Computes the variance of the best split point for two normal distributions.
    The prototype is in the Mathematica file StatisticsSplitPoint2.nb
    which has the computations of the variance using the delta method.
    @param n1    number of samples in first distribution
    @param m1    mean of first distribution
    @param v1    variance of first distribution
    @param n2    number of samples in second distribution
    @param m2    mean of second distribution
    @param v2    variance of second distribution
    @param whichSol for which solution to compute the variance
*/
double UnidimensionalQDAVariance(double n1, double m1, double v1,
                                 double n2, double m2, double v2, int whichSol)
{

    whichSol=6;

#ifdef DEBUG_PRINT

    std::cout << "UnidimentionalQDAVariance " << n1 << " " << n2 << " " << v1 << " " << v2;
    std::cout << " " << m1 << " " << m2 << "\t";
#endif

    double variance;

    if (v1<SMALL_POZ_VALUE  || v2<SMALL_POZ_VALUE )
    {
        if (n1<3 || n2<3)
            variance=pow2(m1-m2);
        else
            variance=SMALL_POZ_VALUE;
        return variance;
    }

    if (std::isnan(m1) || std::isnan(m2) || std::isnan(v1) || std::isnan(v2))
    {
        variance=LARGE_POZ_VALUE;
        return variance;
    }


    if ( IsEqual(v1,v2) && IsEqual(m1,m2) )
        variance = SMALL_POZ_VALUE;
    else
    {
        double CLog=log(n1/n2*sqrt(v2/v1));
        double SR=sqrt(v1*v2*( pow2(m1-m2)+2*(-v1+v2)*CLog ));
        double Q1=v1*(pow2(m1-m2)-v1+v2-2*(v1-2*v2)*CLog);
        double Q2=v2*(pow2(m1-m2)+v1-v2-2*(2*v1-v2)*CLog);
        double Q3=2*(pow(m2,4.0)-2*pow2(m2)*n2*v2+n2*pow2(v2));
        double Q4=2*(pow(m1,4.0)-2*pow2(m1)*n1*v1+n1*pow2(v1));
        double Q5=-m2*v1+m1*v2;

        // we compute all solutions even if not needed. any decent
        // compiler should get rid of unnecessary code

        double Sol1=
            ( v1*pow2(v2) ) / ( n1*pow2(v1-v2) ) * pow2(-1+(m1-m2)*v1/SR)+
            ( v2*pow2(v1) ) / ( n2*pow2(v1-v2) ) * pow2( 1-(m1-m2)*v2/SR)-
            Q4/( n1*(n1-1)*pow(v1-v2,4.0) ) * pow2(  Q5-SR+(v1-v2)*( m2+Q2/SR) )-
            Q3/( n2*(n2-1)*pow(v1-v2,4.0) ) * pow2( -Q5+SR+(v1-v2)*(-m1+Q1/SR) );

        double Sol2=
            ( v1*pow2(v2) ) / ( n1*pow2(v1-v2) ) * pow2( 1+(m1-m2)*v1/SR)+
            ( v2*pow2(v1) ) / ( n2*pow2(v1-v2) ) * pow2(-1-(m1-m2)*v2/SR)-
            Q4/( n1*(n1-1)*pow(v1-v2,4.0) ) * pow2( Q5+SR-(v1-v2)*(-m2+Q2/SR) )-
            Q3/( n2*(n2-1)*pow(v1-v2,4.0) ) * pow2( Q5+SR+(v1-v2)*( m1+Q1/SR) );

        // solution for the first order equation. Here v1=v2
        double QL=v1*log(n1/n2)/pow2(m1-m2);
        double Sol0=2*pow2(m1-m2)*pow2(QL)*(1.0/n1+1.0/n2)+v1/(4*n1)*pow2(1-2*QL)+
                    +v1/(4*n2)*pow2(1+2*QL);

        double SolDefault=(n1*pow2(v1)+n2*pow2(v2))/pow2(n1*sqrt(v1)+n2*sqrt(v2));
        //double SolDefault=(n1*v1+n2*v2)/pow2(n1+n2);

        switch (whichSol) {
        case 0:
            variance=Sol0;
            break;

        case 1:
            variance=Sol1;
            break;

        case 2:
            variance=Sol2;
            break;

        case 3:
            variance=SolDefault;
            break;

        case 5:
            variance=(n1*v1+n2*v2)/pow2(n1+n2);
            break;

        case 6: // paper version
            {
                double s1=sqrt(v1);
                double s2=sqrt(v2);
                variance=s1*s2/(s1+s2)*(s1/n1+s2/n2);
            }
            break;

        default: // very degenerate solution
            variance=LARGE_POZ_VALUE;
        }
#ifdef DEBUG_PRINT
        std::cout << "Sol0=" << Sol0 << "\tSol1=" << Sol1 << "\tSol2=" << Sol2;
        std::cout << "\tSolDefault=" << SolDefault << "\twhichSol=" << whichSol;
#endif

        // assert (isfinite(variance)); //NaN?
    }

#ifdef DEBUG_PRINT
    std::cout << "\tvariance=" << variance << std::endl;
#endif

    return (finite(variance)==true) ? variance : LARGE_POZ_VALUE;
}


/**
   The hyperplane is orthogonal on one of the axis not oblique.
   the best axis to split on is determined by comparing the value of
   the t-test for the directions, i.e.

      (eta1-eta2)^2
   t=-----------------
     sigma1^2+sigma2^2

   the split point is determined by solvind QDA on the unidimentional
   space like for the general case
   @param mass         the total mass of the two distribution
   @param alpha_1    probability to belong to the first distribution
   @param mu1        vector with means for each attribute of points with first class label
   @param S1         covariance mattrix for points with first class label
   @param alpha_2    probability to belong to the second distribution
   @param mu2        vector with means for each attribute of points with second class label
   @param S2         covariance mattrix for points with second class label
   @param SeparatingHyperplane    vector where the best separating hyperplane is returned
   @return best gini gain
*/
double ComputeSeparatingHyperplane_Anova(double mass,
        double alpha_1, arma::vec& mu1, arma::mat& S1,
        double alpha_2, arma::vec& mu2, arma::mat& S2,
        arma::vec& SeparatingHyperplane)
{
    arma::uword maxi;
    double maxt = 0.0;
    int N = mu1.n_elem;

    arma::vec comp = square(mu1 - mu2) / (S1.diag() + S2.diag());

    if (any(comp > 0))
        maxt = comp.max(maxi);
    else
        return 0.0; // none of the directions is any good

    // solve QDA on the projection maxi
    int whichSol; //
    double eta = UnidimensionalQDA(alpha_1, mu1[maxi], S1(maxi+1, maxi+1),
                                   alpha_1, mu2[maxi], S2(maxi+1, maxi+1), whichSol);
    double variance = UnidimensionalQDAVariance(mass * alpha_1, mu1[maxi], S1(maxi, maxi),
                                                mass * alpha_1, mu2[maxi], S2(maxi, maxi), whichSol);

    assert(variance > 0.0);

    // form the separating hyperplane
    // by divideing all components by sqrt(variance) we get the
    // quantities computed with the Separating Hyperplane to be
    // centered on 0 and have the standard distribution of the
    // split point 1.0
    SeparatingHyperplane.set_size(N+1);
    double sgn = (mu1[maxi]>eta) ? 1.0 : -1.0;
    SeparatingHyperplane[0] = -sgn * eta / sqrt(variance);
    SeparatingHyperplane[maxi+1] = sgn / sqrt(variance);

    double p11 = alpha_1 * PValueNormalDistribution(mu1[maxi], sqrt(S1(maxi, maxi)), eta);
    double p1_ = p11 + alpha_2 * PValueNormalDistribution(mu2[maxi], sqrt(S2(maxi, maxi)), eta);
    double gini = BinaryGiniGain(p11, alpha_1, p1_);
    return gini;
}

/**
   Compute the quadratic that separates two distributions
   and take the separating hyperplane to be the tangent to it
   in the intersection point with the line between the centers.

   @param mass         the total mass of the two distribution
   @param alpha_1    probability to belong to the first distribution
   @param mu1        vector with means for each attribute of points with first class label
   @param S1         covariance mattrix for points with first class label
   @param alpha_2    probability to belong to the second distribution
   @param mu2        vector with means for each attribute of points with second class label
   @param S2         covariance mattrix for points with second class label
   @param SeparatingHyperplane    vector where the best separating hyperplane is returned
   @return best gini gain
*/
double ComputeSeparatingHyperplane_QDA(double mass,
                                       double alpha_1, arma::vec& mu1, arma::mat& S1,
                                       double alpha_2, arma::vec& mu2, arma::mat& S2,
                                       arma::vec& SeparatingHyperplane)
{

    // Sice we got here we clearly have 2 clusters. So we should better produce a split at least
    // into the regressor space.

    int N = mu1.n_elem;

    /* Check if one of the matrices S1 or S2 has a 0 on the diagonal. If this is the case
       then the variable corresponding to that position is either extremelly predictive (one
       value of the variable determines exactly one of the clusters) or the variable is useless
       (all the values are identical). In the first case we split on it in the seccond we
       substibute the 0 with same small value in both matrices. */

    for (int i=1; i <= N; i++)
    {
        if (S1(i - 1, i - 1) == 0.0 && S2(i - 1, i - 1) == 0.0)
        {
            S1(i - 1, i - 1) = S2(i - 1, i - 1) = SMALL_POZ_VALUE;
            continue;
        }
        if (S1(i - 1, i - 1) == 0.0 || S2(i - 1, i - 1) == 0.0)
        {
            /* split on variable i exclusivelly. Do the split between the projections of the
               centers of the clusters on this direction proportional with the sizes of the clusters.
            */

            SeparatingHyperplane = arma::zeros<arma::vec>(N + 1);
            SeparatingHyperplane[0] = -(alpha_1 * mu1[i - 1] + alpha_2 * mu2[i - 1]);
            SeparatingHyperplane[i] = 1;

            return 1 - pow2(alpha_1) - pow2(alpha_2);
        }
    }

    // compute the Cholesky of S1 in lower part of S1
    if (!chol(S1,S1) || !chol(S2,S2))
    {
        // one of the matrices S1 and S2 is badly conditioned
        std::cerr << "Badly conditioned matrices in ComputeSeparatingHyperplane" << std::endl;

    }

    // compute the square root of determinant of S1 and S2
    double detS1 = det(S1);
    double detS2 = det(S2);

    /* compute the quadratic that separates the distributions
       and take the separating hyperplane to be the tangent to it
       in the intersection point with the line between the centers */

    // Compute C
    double C = 2.0 * log(detS2 / detS1) + log(alpha_1 / alpha_2);

    arma::vec dif = mu2 - mu1;
    arma::vec G1inv_dif = S1.i() * dif;
    double I1 = arma::dot(G1inv_dif, G1inv_dif);

    arma::vec G2inv_dif = S1.i() * dif;
    double I2 = arma::dot(G2inv_dif, G2inv_dif);

    // compute t as solution of (1-t)^2 I1-t^2 I2=C
    double t, delta;
    if (IsZero(I1 - I2))
    {
        // Liner equation
        t = (I1 - C) / (2 * I1);
    }
    else
    {
        // quadratic equation
        delta= C * (I1-I2) + I1 * I2;
        if (delta < 0.0)
        {
            t = -1.0;
        }
        else
        {
            t = (I1 - sqrt(delta)) / ( I1-I2);
            if (t < 0.0 || t > 1.0)
                t = (I1 + sqrt(delta)) / (I1 - I2);
        }
    }

    arma::vec n; // the normal to the Separating Hyperplane

    /* check if the solution makes any sense. If not take a point proportional
       with the sizes (alpha) on the line that goes through the two centers
    */
    if (t < 0.0 || t > 1.0)
    {
        // fix t
        t = alpha_1 / (alpha_1 + alpha_2); // in case alphas are not normalized
        n = normalise(dif);
    }
    else
    {
        // compute n as the tangent to the the quadratic in point mu
        n = normalise(S1.i() * G1inv_dif * (t - 1) - S2.i() * G2inv_dif * t);
    }

    // compute the point which the  Separating Hyperplane crosses
    arma::vec mu = t * mu1 + (1-t) * mu2;

    // compute etai, vari
    double eta1 = arma::dot(n,mu1);
    double eta2 = arma::dot(n,mu2);

    // compute vari=n^T*S1*n directly into vari
    double var1 = arma::as_scalar(n.t() * S1 * n), var2 = arma::as_scalar(n.t() * S2 * n);

    double variance = UnidimensionalQDAVariance(mass * alpha_1, eta1, var1,
                                                mass * alpha_2, eta2, var2, 5);

    assert(variance > 0.0);

    // compute the equation of the SeparatingHyperplane
    double eta = arma::dot(n,mu);
    double sgn = (eta1 > eta ? 1.0 : -1.0);
    SeparatingHyperplane.set_size(N + 1);
    SeparatingHyperplane[0] = -sgn * eta / sqrt(variance);
    SeparatingHyperplane.subvec(1, N) = sgn / sqrt(variance) * n;

    double p11 = alpha_1 * PValueNormalDistribution(mu1, S1, n, mu);
    double p1_ = p11 + alpha_2 * PValueNormalDistribution(mu2, S2, n, mu);
    double gini = BinaryGiniGain(p11, alpha_1, p1_);
    return gini;
}


/**
   The normal of the hyperplane is the best separating direction, i.e.
   the vector on which the projection of the two gaussians is as separated
   as possible as measured by Fisher's discriminant.

   The equations are:
   - normal: n=(alpha_1*S1+alpha_2*S2)^{-1}(mu1-mu2)
   - separating equations: n^T*x-eta=0
   - eta is solution of the equation

   eta^2(1/var1 - 1/var2)-2eta(eta1/var1-eta2/var2)+eta1^2/var1-eta2^2/var2
   = 2ln(alpha1/alpha2)-ln(var1/var2)

   with etai=n^T*mu1 and vari=n^T*Si*n

   S1 and S2 are modified and have their respective Cholesky factorization in them

   @param mass         the total mass of the two distribution
   @param alpha_1    probability to belong to the first distribution
   @param mu1        vector with means for each attribute of points with first class label
   @param S1         covariance mattrix for points with first class label
   @param alpha_2    probability to belong to the second distribution
   @param mu2        vector with means for each attribute of points with second class label
   @param S2         covariance mattrix for points with second class label
   @param SeparatingHyperplane    vector where the best separating hyperplane is returned
   @return best gini gain
*/
double ComputeSeparatingHyperplane_LDA(double mass,
                                       double alpha_1, arma::vec& mu1, arma::mat& S1,
                                       double alpha_2, arma::vec& mu2, arma::mat& S2,
                                       arma::vec& SeparatingHyperplane)
{

    if (alpha_1 < TNNearlyZero || alpha_2 < TNNearlyZero)
        return 0.0; // there's nothing we can do. One of the clusters is nonexistent

    int N = S1.n_rows;
    arma::vec n = mu1 - mu2; // normal to the hyperplane
    double eta1, eta2;
    double var1, var2;

    arma::mat Sw = alpha_1 * S1 + alpha_2 * S2;

    // Solve Sw^{-1}*n into n
    // First compute Cholesky of Sw into lower Sw
    if (!chol(Sw, Sw))
    {
        std::cerr << "Badly conditioned matrices in ComputeSeparatingHyperplane" << std::endl;
        std::cerr << "We use normal anova analysis to find the Separating Hyperplane" << std::endl;
        return ComputeSeparatingHyperplane_Anova(mass,alpha_1, mu1, S1, alpha_2, mu2, S2, SeparatingHyperplane);
    }

    n = normalise(trimatu(Sw.i() * n));
    if (arma::all(n == 0))
        return 0.0; // nothing we can do. Clusters are identical

    // compute etai, vari
    eta1 = arma::dot(n, mu1);
    eta2 = arma::dot(n, mu2);

    // compute vari=n^T*S1*n directly into vari
    var1 = arma::as_scalar(n.t() * S1 * n);
    var2 = arma::as_scalar(n.t() * S2 * n);
    var1 = 0.0;
    var2 = 0.0;

    int whichSol;
    double eta = UnidimensionalQDA(alpha_1, eta1, var1, alpha_2, eta2, var2, whichSol);
    double variance = UnidimensionalQDAVariance(mass * alpha_1, eta1, var1,
                                                mass * alpha_2, eta2, var2, whichSol);

    assert(variance > 0.0);

    // form the separating hyperplane
    // by divideing all components by sqrt(variance) we get the
    // quantities computed with the Separating Hyperplane to be
    // centered on 0 and have the standard distribution of the
    // split point 1.0
    double sgn = (eta1 > eta ? 1.0 : -1.0);
    SeparatingHyperplane.set_size(N + 1);
    SeparatingHyperplane[0] = -sgn * eta / sqrt(variance);
    SeparatingHyperplane.subvec(1, N) = sgn / sqrt(variance) * n;

    arma::vec mu(N);
    mu = eta * n;


    // compute the Cholesky of S1 in lower part of S1
    if (!arma::chol(S1,S1) || !arma::chol(S2,S2))
    {
        // one of the matrices S1 and S2 is badly conditioned
        std::cerr << "Badly conditioned matrices in ComputeSeparatingHyperplane" << std::endl;
        return 0.0; // worst gini
    }

    double p11 = alpha_1 * PValueNormalDistribution(mu1, S1, n, mu);
    double p1_ = p11 + alpha_2 * PValueNormalDistribution(mu2, S2, n, mu);
    double gini = BinaryGiniGain(p11, alpha_1, p1_);
    return gini;
}


}

#endif //_CLUS_SPLITPOINTCOMPUTATION_H_
