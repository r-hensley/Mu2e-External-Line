//
//
//  This routine reads the DST NTuple which was created in a particular
//  g4bl script, and generates a subsequent input file with only the surviving
//  particles
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

#define DST_cxx
#include <DST.h>

//Dummy definition
void DST::Loop() {
}

using namespace std;

//
// If emittance is negative, then it generates a uniform emittance.  Otherwise, it
// generates a Gaussian emittance
//

void Survive(string inputFileName, string outputFileName, 
   int startLoc=0, int nEvents=0) {
// determine the operating mode



// Open the Ntuple and load the DST NTuple
TFile *inputFile = new TFile(inputFileName.c_str());
TTree *t = (TTree*) inputFile->FindObjectAny("DST");
DST *d = new DST(t);

TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
outputFile->mkdir("VirtualDetector");
outputFile->cd("VirtualDetector");
TNtuple *outputNtuple = new TNtuple("tracklist","Track List","x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight");
  


// Loop over the entries in the DST
Long64_t nentries = d->fChain->GetEntriesFast();

int nTotal=0;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if((nEvents>0)&&(jentry>nEvents)) break;

      nTotal++;
      Long64_t ientry = d->LoadTree(jentry);
      if (ientry < 0) break;
      d->fChain->GetEntry(jentry);  
      // Calculate the emittance at the point of generation
      // If we're here, then it has passed the emittance cut so we load it into
      // the output Ntuple
      Float_t x,y,z,Px,Py,Pz,PDGid,EventID,Weight;
      // Choose starting location based on startLoc
      if(startLoc==0) {            // beginning of beam line 
		  x = d->Z2201_x;
		  y = d->Z2201_y;
		  z = d->Z2201_z;
		  Px = d->Z2201_Px;
		  Py = d->Z2201_Py;
		  Pz = d->Z2201_Pz;
		  PDGid = d->Z2201_PDGid;
		  EventID = d->Z2201_EventID;
		  Weight = d->Z2201_Weight;
      }else if (startLoc==2) {     // Middle of AC Dipole
		  x = d->Z137004_x;
		  y = d->Z137004_y;
		  z = d->Z137004_z;
		  Px = d->Z137004_Px;
		  Py = d->Z137004_Py;
		  Pz = d->Z137004_Pz;
		  PDGid = d->Z137004_PDGid;
		  EventID = d->Z137004_EventID;
		  Weight = d->Z137004_Weight;
      } else if (startLoc==3) {    // After final collimator
		  x = d->Z160749_x;
		  y = d->Z160749_y;
		  z = d->Z160749_z;
		  Px = d->Z160749_Px;
		  Py = d->Z160749_Py;
		  Pz = d->Z160749_Pz;
		  PDGid = d->Z160749_PDGid;
		  EventID = d->Z160749_EventID;
		  Weight = d->Z160749_Weight;
      }  else {
          cerr << "Unknown start location code "<<startLoc<<endl;
          break;
      }      
      outputNtuple->Fill(x,y,z,Px,Py,Pz,0.,PDGid,EventID,0.,0.,Weight); 
   }
   
cout << "Processed "<<nTotal<<" events."<<endl;

outputFile->Write();
outputFile->Close();

}
