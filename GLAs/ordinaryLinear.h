#include <armadillo>
#include <math.h>
#include <boost/math/distributions.hpp>
#include "arma/Types/VECTOR.h"
#include "arma/Types/MATRIX.h"
#include "base/Types/INT.h"
#include "base/Types/BIGINT.h"

#define ADD(x) x += other.x;

/*
 *  This is an implementation of Ordinary Linear Regression
 */

/** Information for Meta-GLAs
 *  GLA_DESC
 *
 *  NAME(</OrdinaryLinear/>)
 *  INPUTS(</(x, VECTOR), (y, DOUBLE)/>)
 *  OUTPUTS(</(count, BIGINT)/>)
 *  CONSTRUCTOR(</(width, BIGINT)/>)
 *  RESULT_TYPE(</single/>)
 *
 *  LIBS(armadillo)
 *
 *  END_DESC
 *
 */

using namespace arma;

// Declaration
class OrdinaryLinear;

class OrdinaryLinear {

    uint64_t numRows;           // number of rows processed
    MATRIX sparsity;            // X^T times X
    VECTOR rhs;                 // x times y, right hand side of equation
    VECTOR beta;                // vector or parameters for estimatio

public:

    OrdinaryLinear( const INT  & width )
				:	numRows         (0),
					rhs             (width),
					sparsity        (width, width),
					beta            (width) {
      beta.zeros();
      rhs.zeros();
      sparsity.zeros();
    }

    void AddItem( const VECTOR & x, const DOUBLE & y ) {
      ++numRows;
      sparsity += x * trans(x);
			rhs += x * y;
    }

    void AddState( const OrdinaryLinear  & other ) {
      ADD(sparsity);
      ADD(rhs);
      ADD(numRows);
    }

    void CalcBeta() {
      solve(beta, adjSparsity, adjRHS);
    }

    void PrintBeta() {
      cout << endl;
      cout << beta << endl;
    }

    void GetResult( BIGINT & count ) {
			CalcBeta();
      count = 1;
      PrintBeta();
    }
};
