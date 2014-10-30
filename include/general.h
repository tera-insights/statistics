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

#if !defined _GENERAL_H
#define _GENERAL_H

// uncoment the following if debugin
#define    DEBUGON

// uncomment to print lots of debug messages
// #define DEBUG_PRINT

// uncomment the following if XML is used
// #define CLUS_USE_XML

//#define TNT_BOUNDS_CHECK

#if defined(DEBUGON)
#include <iostream>
#endif

// constants used throught the program
#define MaxIterations        50
#define MAXREAL              1.0e+100
#define MINREAL              -1.0e+100
#define BLKSIZE              4096
#define BIG_NR               1000000
#define SMALL_POZ_VALUE      3.0E-13
#define LARGE_POZ_VALUE      3.0E+100
#define MAX_CATEGORIES       10000
#define MAX_VARIABLES        1000

#if !defined INT_MAX
#define INT_MAX              2147483647
#endif

/* error values return by some functions */
#define    EBADDIM            1

// robost functions to check equality and closedness to zero
#define NonEqual(X,Y)         ( (fabs(X)+fabs(Y))==0.0?0:NonZero((X-Y)/(fabs(X)+fabs(Y))) )
#define IsEqual(X,Y)          ( (fabs(X)+fabs(Y))==0.0?1:IsZero((X-Y)/(fabs(X)+fabs(Y))) )
#define NonZero(X)            ( ( (X)>TNNearlyZero ) || ( (X)<-TNNearlyZero ) )
#define IsZero(X)             ( ( (X)<TNNearlyZero ) && ( (X)>-TNNearlyZero ) )
#define NonZeroM(X)           ( ( (X)>TNNearlyZero*100 ) || ( (X)<-TNNearlyZero*100 ) )

#define RANDOM01FLOAT         ( (1.0*rand())/RAND_MAX )

#if defined(DEBUGON)
#define ExitIf( Cond, Text )     if ( Cond ) { \
                                    std::cout <<  Text << std::endl; \
                                    exit(1); \
                                 }
#else
#define ExitIf( Cond, Text )
#endif

/// Computes second power of a double. Very useful
inline double pow2(double x)
{
    return x*x;
}

#endif /* _GENERAL_H */
