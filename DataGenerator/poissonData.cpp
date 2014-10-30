#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/datapath/statistics/Data/poissonData.txt");
  random_device rd;
  mt19937 rand (rd());
  double intercept = generate_canonical<double, 40>(rand);
  double b1 = generate_canonical<double, 40>(rand);
  double b2 = generate_canonical<double, 40>(rand);
  double b3 = generate_canonical<double, 40>(rand);
  myfile << intercept << "," << b1 << "," << b2 << "," << b3 << endl;
  for(int count = 0; count < 1000000; count ++){
    double x1 = generate_canonical<double, 40>(rand);
    double x2 = generate_canonical<double, 40>(rand);
    double x3 = generate_canonical<double, 40>(rand);
    double param = std::exp(intercept + x1 * b1 + x2 * b2 + x3 * b3);
    std::poisson_distribution<int> distribution (param);
    double randomVar = distribution(rand);
    /*double p = generate_canonical<double, 40>(rand);
    int counter;
    int currentFactorial = 1;
    for(counter = 0; p > 0; counter ++){
      double currentProbability = pow(param, counter) * exp(-param) / currentFactorial;
      p -= currentProbability;
      currentFactorial *= (counter + 1);
    }*/
    myfile << x1 << "," << x2 << "," <<  x3 << "," << randomVar << endl;
  } 
  cout << intercept << " " << b1 << " " << b2 << " " << b3 << endl;
  myfile.close();
  return 0;
}
