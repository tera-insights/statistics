#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/datapath/statistics/Data/seqData.txt");
  for(int counter = 0; counter < 100; counter ++){
    myfile << (counter + 1) << endl;
  }
  myfile.close();
  return 0;
}
