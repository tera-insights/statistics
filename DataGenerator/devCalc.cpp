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
  double dev = 0; //deviance
  int n = 10000; //count
  int p = 2; //family specific, 0 for gaussian
  double x1, x2, x3, y, u;
  for(int counter = 0; counter < n; counter ++)
    {
      getline(myfile, s, ',');
      x1 = atof(s.c_str());
      getline(myfile, s, ',');
      x2 = atof(s.c_str());
      getline(myfile, s, ',');
      x3 = atof(s.c_str());
      getline(myfile, s, ',');
      y = atof(s.c_str());
      getline(myfile, s);
      u = 1/(intercept);
      dev += 2 * (std::log(u / y) + (y - u)/u); // sum (y-u)^2 / V(u)
    }
  cout << dev
;
}
