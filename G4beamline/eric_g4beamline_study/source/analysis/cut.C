//
// Very simple program to make a momentum cut on a g4Beamline output file
//
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>

#define G4BLVirtualDetector_cxx
#include "G4BLVirtualDetector.h";

using namespace std;


Float_t rProj(Float_t zProj, Float_t x, Float_t y, Float_t z, Float_t Px, Float_t Py, Float_t Pz);
Float_t quadAdd(Float_t x, Float_t y);

TTree *inTree=NULL;
TFile *inputFile=NULL;
TFile *outputFile=NULL;

TNtuple *outputNtuple=NULL;
TNtuple *analysisNtuple=NULL;

Float_t x,y,z,Px,Py,Pz,r,pT,pTot,angle;

G4BLVirtualDetector *d=NULL;



TDatabasePDG pdb;

void cut(string fn,float pCut) {

string inputFileName=fn;
string outputFileName="cut-"+fn;

cout << "Input File Name:" << inputFileName <<endl;
cout << "Output File Name:" << outputFileName <<endl;
cout << "Momentum Cut:" << pCut <<endl;

outputFile = new TFile(outputFileName.c_str(),"RECREATE");
// Output Ntuple for analysis
analysisNtuple = new TNtuple("anaTracks","Analysis Tracks","x:y:z:Px:Py:Pz:r:pT:pTot:angle:dFlag:PDGid:beta:q");

// Input Ntuple for next part of analysis
outputFile->mkdir("VirtualDetector");
outputFile->cd("VirtualDetector");
outputNtuple = new TNtuple("cutTracks","Cut Track List","x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight");

TFile f(inputFileName.c_str());
inTree = (TTree *) f.FindObjectAny("DetC1B");
d = new G4BLVirtualDetector(inTree);


// Loop over input tree and cut low momentum tracks

Long64_t nentries = d->fChain->GetEntriesFast(); 
Long64_t nbytes = 0, nb = 0;
float dFlag=0.;  //Detector flag.  Not currently used.
for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = d->LoadTree(jentry);
    if (ientry < 0) break;

    int count = jentry+1;
    int place = (int) log10(count);
    int scale = power(10,place);
    if((count%scale)==count) count << "Processing Track "<<count<<" out of "<<nentries<<endl;  

    nb = d->fChain->GetEntry(jentry);   nbytes += nb;
    x=d->x;
    y=d->y;
    z=d->z;
    Px=d->Px;
    Py=d->Py;
    Pz=d->Pz;
    r = quadAdd(x,y);
    pT = quadAdd(Px,Py);
    pTot = quadAdd(pT,Pz);
    if(pTot<pCut) continue; //Cut low momentum tracks
    
    angle = asin(pT/pTot);
    float beta=0;
    int q=0;
    TParticlePDG *p = pdb.GetParticle((int) d->PDGid);
    if(p!= NULL) {
       float M = 1000.*p->Mass();
       beta=pTot/quadAdd(pTot,M);
       q = p->Charge();
    }
    // Fill the Ntuples    
    analysisNtuple->Fill(x,y,z,Px,Py,Pz,r,pT,pTot,angle,dFlag,d->PDGid,beta,(Float_t) q);  //Store to output NTuple
    outputNtuple->Fill(x,y,z,Px,Py,Pz,d->t,
     d->PDGid,d->EventID,d->TrackID,d->ParentID,d->Weight);

}
//
// Close output files
//
outputFile->Write();
outputFile->Close();

}
Float_t quadAdd(Float_t x, Float_t y) {
   return(sqrt(x*x+y*y));
}

