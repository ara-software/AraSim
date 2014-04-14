#include "Efficiencies.h"
#include <string>
using namespace std;
Efficiencies::Efficiencies(int nRx,string outputdir) {
  string stemp=outputdir+"/output.txt";
  fout.open(stemp.c_str());
  l1Counter.push_back(nRx);
}


Efficiencies::~Efficiencies() {

    l1Counter.clear();

}


void Efficiencies::summarize() {

  fout << "Summarize results.\n";
  fout.close();


}
void Efficiencies::incrementL1Counter(vector<int> l1Hits) {
  for (int i=0;i<(int)l1Hits.size();i++) {
    l1Counter.push_back(l1Hits[i]);
  }

}
