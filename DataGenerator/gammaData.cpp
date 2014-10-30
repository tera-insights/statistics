#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/datapath/statistics/Data/gammaData.txt");
  random_device rd;
  mt19937 rand (rd());
  double intercept = generate_canonical<double, 40>(rand) + 1;
  double b1 = generate_canonical<double, 40>(rand) + 1;
  double b2 = generate_canonical<double, 40>(rand) + 1;
  double b3 = generate_canonical<double, 40>(rand) + 1;
  for(int counter = 0; counter < 10000; counter ++){
    //double shape = generate_canonical<double, 40>(rand) + rand() % 9 + 1;
    int shape = 5;
    double x1 = generate_canonical<double, 40>(rand) + 1;
    double x2 = generate_canonical<double, 40>(rand) + 1;
    double x3 = generate_canonical<double, 40>(rand) + 1;
    double mu = 1/(b1 * x1 + b2 * x2 + b3 * x3 + intercept);
    double scale = mu / shape;
    std::gamma_distribution<double> distribution(shape, scale);
    double randVar = distribution(rand);
    myfile << x1 << "," << x2 << "," <<  x3 << "," << randVar << endl;
  }
  cout << intercept << " " << b1 << " " << b2 << " " << b3 << endl;
  myfile.close();
  return 0;
}
