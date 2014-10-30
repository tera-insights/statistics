#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/datapath/statistics/Data/factorData.txt");
  random_device rd;
  mt19937 rand (rd());
  double intercept = generate_canonical<double, 40>(rand) + 1;
  double b1 = generate_canonical<double, 40>(rand) + 1;
  double b2 = generate_canonical<double, 40>(rand) + 1;
  double b3 = generate_canonical<double, 40>(rand) + 1;
  double a1 = generate_canonical<double, 40>(rand) + 1;
  double a2 = generate_canonical<double, 40>(rand) + 1;
  double a3 = generate_canonical<double, 40>(rand) + 1;
  for(int counter = 0; counter < 10000; counter ++){
    //double shape = generate_canonical<double, 40>(rand) + rand() % 9 + 1;
    int shape = 5;
    double x2 = generate_canonical<double, 40>(rand);
    double x1 = generate_canonical<double, 40>(rand);
    double factorInput;
		double factorInput2;
    char factorType;
		char factorType2;
    if (x1 < .33) {
      factorInput = b1;
      factorType = 'a';
    } else if (x1 < .66) {
      factorInput = b2;
      factorType = 'b';
    } else {
      factorInput = b3;
      factorType = 'c';
    }
    if (x2 < .33) {
      factorInput2 = a1;
      factorType2 = 'a';
    } else if (x2 < .66) {
      factorInput2 = a2;
      factorType2 = 'b';
    } else {
      factorInput2 = a3;
      factorType2 = 'c';
    }
    double mu = 1/(factorInput + factorInput2 + intercept);
    double scale = mu / shape;
    std::gamma_distribution<double> distribution(shape, scale);
    double randVar = distribution(rand);
    myfile << factorType << "," << factorType2 << "," << randVar << endl;
  }
  cout << intercept << " " << b1 << " " << b2 << " " << b3 << " " << a1 << " " << a2 << " " << a3 << endl;
  myfile.close();
  return 0;
}
