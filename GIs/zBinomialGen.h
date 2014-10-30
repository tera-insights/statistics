// To be included by module
#include <sstream>
#include <vector>
#include <string>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include "GIStreamInfo.h"

/*
 * GI_DESC
 *  NAME(BinomialGen)
 *  OUTPUTS(</(input1, DOUBLE), (input2, DOUBLE), (input3, DOUBLE), (y1, INT), (y2, INT)/>)
 *
 * END_DESC
 */

class BinomialGen {

    // Template parameter
    int counter;
    int width;
    double intercept;
    std::mt19937 rand;
    VECTOR  xCoef;
    int numTrial;

 public:

 BinomialGen(GIStreamProxy& stream) :
    counter	(0),
      width	(3),
      rand	(stream.get_id()),
      xCoef	(3),
      intercept	(generate_canonical<double, 40>(rand)),
      numTrial	(rand() % 100)
      {
	cout << "intercept " << intercept << endl;
	for(int counter = 0; counter < width; counter++){
	  xCoef(counter) = generate_canonical<double,40>(rand);
	  cout << "coefficient " << counter << ": " << xCoef(counter) << endl;
	}
	
      }

    bool ProduceTuple(double& x1, double& x2, double& x3, int& y1, int& y2) {
        if( counter < 1000000 ) {
	  std::ofstream myFile;
	  myFile.open("/home/jon//datapath/statistics/Data/binomData.txt", std::ios_base::out | std::ios_base::app);
	  x1 = generate_canonical<double, 40>(rand);
	  x2 = generate_canonical<double, 40>(rand);
	  x3 = generate_canonical<double, 40>(rand);
	  double probability = 1/(1 + std::exp(-(intercept + x1 * xCoef(0) + x2 * xCoef(1) + x3 * xCoef(2))));
	  /*
	  x.setSize(width);
	  for(int counter = 0; counter < width; counter ++){
	    x(counter) = generate_canonical<double, 40>(rand);
	  }
	  double probability = intercept + dot(x, xCoef);
	  */
	  std::binomial_distribution<int> binom(numTrial, probability);
	  y1 = binom(rand);
	  y2 = numTrial - y1;
	  ++counter;
	  std::ostringstream ss;
	  ss << x1 << "," << x2 << "," << x3  << "," << y1 << "," << y2 << std::endl;
	  //std::string output = x1 + ", " + x2 + ", " + x3 + ", " + y1 + ", " + y2 + "\n";
	  myFile << ss.str();
	  myFile.close();

	  return true;
        } 
	else 
	  return false;
    }
};
