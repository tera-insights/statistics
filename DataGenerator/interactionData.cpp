#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/datapath/statistics/Data/interactionData.txt");
  random_device rd;
  mt19937 rand (rd());
  double intercept = generate_canonical<double, 40>(rand) + 1;
  double beta [9];
  for(int counter = 0; counter < 9; counter ++)
    beta[counter] = generate_canonical<double, 40>(rand) + 1;
  for(int counter = 0; counter < 10000; counter ++){
    //double shape = generate_canonical<double, 40>(rand) + rand() % 9 + 1;
    int shape = 5;
    int x2 = rand() % 3;
    int x1 = rand() % 3;
    double factorInput;
    char factorType1;
    char factorType2;
    switch (x2) {
      case 0:
        factorType2 = 'a';
        break;
      case 1:
        factorType2 = 'b';
        break;
      case 2:
        factorType2 = 'c';
        break;
    }
    switch (x1) {
      case 0:
        factorType1 = 'a';
        break;
      case 1:
        factorType1 = 'b';
        break;
      case 2:
        factorType1 = 'c';
        break;
    }
    factorInput = beta[3 * x1 + x2];
    double mu = 1/(factorInput);
    double scale = mu / shape;
    std::gamma_distribution<double> distribution(shape, scale);
    double randVar = distribution(rand);
    myfile << factorType1 << "," << factorType2 << "," << randVar << endl;
  }
  cout << "beta data: " << endl;
  for(int counter = 0; counter < 9; counter ++)
    cout << beta[counter] << endl;
  myfile.close();
  return 0;
}
