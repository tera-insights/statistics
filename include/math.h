/*

Copyright (c) 2014, Tera Insights, LLC
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

#ifndef _CLUS_MATH_H_
#define _CLUS_MATH_H_

#include <limits>
#include <cmath>

namespace CLUS {

// Math constants used throughout the programs
constexpr const double LAMBDA_PREC = 1.0E-5;
constexpr const double PREC_ARITH = 1.0E-7;
constexpr const double TNNearlyZero = 1.0E-5;
constexpr const double EPSILON = 1.0E-5;
constexpr const double JacobiTolerance = 1.0E-6;

// Useful functions for comparing floating point numbers

/**
 * Determine if x and y are equal with a given epsilon.
 *
 * x and y are considered equal if the absolute difference between the two is
 * at most epsilon.
 */
inline
bool is_equal_eps(double x, double y, double eps = EPSILON) {
	return std::fabs(x - y) <= eps;
}

/**
 * Determine if x and y are not equal with a given epsilon.
 *
 * x and y are considered not equal if the absolute difference between the two
 * is more than epsilon.
 */
inline
bool is_not_equal_eps(double x, double y, double eps = EPSILON) {
	return std::fabs(x - y) > eps;
}

/**
 * Determine if x is greater than or equal to y with a given epsilon.
 *
 * x >= y is considered true if (x - y) is at least negative epsilon.
 */
inline
bool is_greater_equal_eps(double x, double y, double eps = EPSILON) {
	return (x - y) >= -eps;
}

/**
 * Determine if x is strictly greater than y with a given epsilon.
 *
 * x > y is considered true if (x - y) is at least epsilon. That is,
 * x and y must differ by at least epsilon, and x must be greater than
 * y.
 */
inline
bool is_greater_eps(double x, double y, double eps = EPSILON) {
	return (x - y) > eps;
}

/**
 * Determine if x is less than or equal to y with a given epsilon.
 *
 * x <= y is considered true if (x - y) is at most positive epsilon.
 */
inline
bool is_less_equal_eps(double x, double y, double eps = EPSILON) {
	return (x - y) <= eps;
}

/**
 * Determine if x is strictly less than y with a given epsilon.
 *
 * x < y is considered true if (x - y) is at most negative epsilon.
 */
inline
bool is_less_eps(double x, double y, double eps = EPSILON) {
	return (x - y) < -eps;
}

}

#endif // _CLUS_MATH_H_
