#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main(){
  ifstream myfile ("/home/jon/datapath/statistics/Data/gammaData.txt");
  double intercept, b1, b2, b3;
  string s;
  getline(myfile, s, ' ');
  intercept = atof(s.c_str());
  getline(myfile, s, ' ');
  b1 = atof(s.c_str());
  getline(myfile, s, ' ');
  b2 = atof(s.c_str());
  getline(myfile, s);
  b3 = atof(s.c_str());
  double phi = 0; //disp param
  int n = 10000; //count
  int p = 2; //family specific, 0 for gaussian
  double x1, x2, x3, y1;
  for(int counter = 0; counter < n; counter ++)
    {
      getline(myfile, s, ',');
      x1 = atof(s.c_str());
      getline(myfile, s, ',');
      x2 = atof(s.c_str());
      getline(myfile, s, ',');
      x3 = atof(s.c_str());
      getline(myfile, s);
      y1 = atof(s.c_str());
      double mu = 1/(intercept + x1 * b1 + x2 * b2 + x3 * b3);
      phi += std::pow(y1 - mu, 2) / std::pow(mu, p); // sum (y-u)^2 / V(u)
    }
  phi /= n - p;
  cout << phi;
}
