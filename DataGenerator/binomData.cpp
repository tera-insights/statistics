#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/datapath/statistics/Data/binomData.txt");
  random_device rd;
  mt19937 rand (rd());
  double intercept = generate_canonical<double, 40>(rand);
  double b1 = generate_canonical<double, 40>(rand);
  double b2 = generate_canonical<double, 40>(rand);
  double b3 = generate_canonical<double, 40>(rand);
  myfile << intercept << " " << b1 << " " <<  b2 << " " << b3 << endl;
  for(int counter = 0; counter < 10000; counter ++){
    int numTrial = rand() % 100 + 1;
    double x1 = generate_canonical<double, 40>(rand);
    double x2 = generate_canonical<double, 40>(rand);
    double x3 = generate_canonical<double, 40>(rand);
    double alpha = 1/(1 + std::exp(-(intercept + x1 * b1 + x2 * b2 + x3 * b3)));
    std::binomial_distribution<int> binom(numTrial, alpha);
    int y1 = binom(rand);
    int y2 = numTrial - y1;
    myfile << x1 << "," << x2 << "," <<  x3 << "," << y1 << "," << y2 <<endl;
    }
  cout << intercept << " " << b1 << " " << b2 << " " << b3 << endl;
  myfile.close();
  return 0;
}
