#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main () {
  ofstream myfile ("/home/jon/Data/randIntData.txt");
  for(int counter = 0; counter < 100000; counter ++){
    myfile << (double) (rand() % 10000) / 10000 << ","
        << (double) (rand() % 10000) / 10000 << ","
        << (double) (rand() % 10000) / 10000 << ","
        << (double) (rand() % 10000) / 10000 << "," << endl;
  }
  myfile.close();
  return 0;
}
