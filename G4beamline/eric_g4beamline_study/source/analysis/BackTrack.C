//
//
// This program will take the tracks that make it to the target in g4beamline and backtrack them
//
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <TNtuple.h>

#include <math.h>

#define DST_cxx
#include "DST.h"

using namespace std;
Float_t quadAdd(Float_t x, Float_t y);

DST *dst=NULL;

Float_t rCut=3.15;

void BackTrack(string fileName) {

Float_t x,y,r;


cout << "Input file name: "<<fileName<<endl;
string outputFileName = "backtrack-"+fileName;
cout << "Output file name: "<<outputFileName<<endl;
string eventFileName = "events-"+fileName+".txt";
cout << "Event file list: "<<eventFileName<<endl;



TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
outputFile->mkdir("VirtualDetector");
outputFile->cd("VirtualDetector");
TNtuple *outputNtuple = new TNtuple("backtrack","Backward Track List","x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight");
  

TFile f(fileName.c_str());

ofstream eventFile;
eventFile.open(eventFileName.c_str());

TTree *dstTree = (TTree *) f.FindObjectAny("DST");
dst = new DST(dstTree);


Long64_t nentries = dst->fChain->GetEntriesFast(); 
Long64_t nbytes = 0, nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = dst->LoadTree(jentry);
    if (ientry < 0) break;
    nb = dst->fChain->GetEntry(jentry);   nbytes += nb;
// Load the position at the target
    x=dst->Z244786_x;
    y=dst->Z244786_y;
   r = quadAdd(x,y);
// Write it out if it hit the target
    if(r<rCut)  { outputNtuple->Fill(dst->Z135163_x,dst->Z135163_y,dst->Z135163_z,
     dst->Z135163_Px,dst->Z135163_Py,dst->Z135163_Pz,dst->Z135163_t,
     dst->Z135163_PDGid,dst->Z135163_EventID,dst->Z135163_TrackID,dst->Z135163_ParentID,dst->Z135163_Weight);
     int EventID= (int) dst->Z135163_EventID;
     eventFile<< EventID << endl;
     }
}
//
// Close output files
//
outputFile->Write();
outputFile->Close();
eventFile.close();

}
Float_t quadAdd(Float_t x, Float_t y) {
   return(sqrt(x*x+y*y));
}
