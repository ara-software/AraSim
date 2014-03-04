////////////////////////////////////////////////////////////////////////////////////////////////
//class Efficiencies:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef EFFICIENCIES_H
#define EFFICIENCIES_H

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

using namespace std;
class Efficiencies {


public:
  
  Efficiencies(int nRx,string outputdir);
  ~Efficiencies();

  void incrementL1Counter(vector<int> l1Hits);
  ofstream fout;
  void summarize();
  vector<int> l1Counter;

protected:


};

#endif //EFFICIENCIES_H
