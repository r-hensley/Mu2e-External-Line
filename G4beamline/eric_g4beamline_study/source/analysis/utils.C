//
// Some handy utility programs
//
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <TNtuple.h>
#include <TRandom2.h>
#include <TFile.h>

#include <math.h>

TNtuple *getNtuple(TFile *,string);
int getEntries(string,string);

TNtuple *getNtuple(TFile *f,string name) {
   return((TNtuple *) f->FindObjectAny(name.c_str()));
}

int getEntries(string fileName, string NtupleName) {
  TFile *f = new TFile(fileName.c_str());
  TNtuple *n = getNtuple(f, NtupleName);
  int nEnt = n->GetEntries();
  cout << "Ntuple " << NtupleName << " in file " << fileName << " has " <<nEnt<<" entries."<<endl;
  delete f;
  return(nEnt);

}