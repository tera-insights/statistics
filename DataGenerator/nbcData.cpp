#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/datapath/statistics/Data/nbcData.txt");
  double x1 = 2.5;
  double y1 = 2.5;
  double x2 = 19.10;
  double y2 = 25.10;
  double x3 = -10.0;
  double y3 = 15.0;
  double sx1 = 2.0;
  double sy1 = 1.5;
  double sx2 = 1.0;
  double sy2 = 1.5;
  double sx3 = 0.5;
  double sy3 = 0.8;
 std::default_random_engine generator;
  std::normal_distribution<double> distribution_x1(x1, sx1);
  std::normal_distribution<double> distribution_y1(y1, sy1);
  std::normal_distribution<double> distribution_x2(x2, sx2);
  std::normal_distribution<double> distribution_y2(y2, sy2);
  std::normal_distribution<double> distribution_x3(x3, sx3);
  std::normal_distribution<double> distribution_y3(y3, sy3);
  for(int counter = 0; counter < 10000; counter ++){
    x1 = distribution_x1(generator);
    y1 = distribution_y1(generator);
    x2 = distribution_x2(generator);
    y2 = distribution_y2(generator);
    x3 = distribution_x3(generator);
    y3 = distribution_y3(generator);
    myfile << "a," << x1 << "," << y1 << endl;
    myfile << "b," << x2 << "," << y2 << endl;
    myfile << "c," << x3 << "," << y3 << endl;
  }
  myfile.close();
  return 0;
}
