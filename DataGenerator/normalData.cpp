#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/datapath/statistics/Data/normalData.txt");
  random_device rd;
  mt19937 rand (rd());
  double intercept = generate_canonical<double, 40>(rand);
  double b1 = generate_canonical<double, 40>(rand);
  double b2 = generate_canonical<double, 40>(rand);
  double b3 = generate_canonical<double, 40>(rand);
  myfile << intercept << " " << b1 << " " << b2 << " " << b3 << endl;
  for(int counter = 0; counter < 10000; counter ++){
    //double shape = generate_canonical<double, 40>(rand) + rand() % 9 + 1;
    int var = 1;
    double x1 = generate_canonical<double, 40>(rand);
    double x2 = generate_canonical<double, 40>(rand);
    double x3 = generate_canonical<double, 40>(rand);
    double mu = (b1 * x1 + b2 * x2 + b3 * x3 + intercept);
    std::normal_distribution<double> distribution(mu, var);
    double randVar = distribution(rand);
    myfile << x1 << "," << x2 << "," <<  x3 << "," << randVar << endl;
  }
  cout << intercept << " " << b1 << " " << b2 << " " << b3 << endl;
  myfile.close();
  return 0;
}
