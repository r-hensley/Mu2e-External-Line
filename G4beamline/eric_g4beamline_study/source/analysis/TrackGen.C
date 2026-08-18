//
//
// This program will generate a track distribution which is uniform in X phase space
// an gaussian in Y phase space
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


using namespace std;
Float_t quadAdd(Float_t x, Float_t y);

// Some important kinematic cuts used for calcul

  double mp=938.272;
  double K=8000;
  double E=K+mp;
  double pMean = sqrt(E*E-mp*mp);
  double betagamma = pMean/mp;
  double betaXAC = 250.176;
  double betaYAC = 4.503;
  double sigmaX;
  double sigmaPx;
  double sigmaY;
  double sigmaPy;


  double PxMean;

// If emittance is negative, then it generates a uniform emittance.  Otherwise, it
// generates a 95% Gaussian emittance

void TrackGen(int nEvents, double emitX=-30.,double emitY=+15.,double zStart=99668.45) {
bool XUniform=false;
bool YUniform=false;

if(emitX<0) XUniform = true;
if(emitY<0) YUniform = true;

emitX=fabs(emitX);
emitY=fabs(emitY);


//Calculate some stuff
// Calculate some kinematic paramters
  sigmaX = sqrt(betaXAC*emitX/(betagamma));
  sigmaY = sqrt(betaYAC*emitY/(betagamma));
  sigmaPx = pMean*sqrt(emitX/(betaXAC*betagamma))/1000;
  sigmaPy = pMean*sqrt(emitY/(betaYAC*betagamma))/1000;
if(!XUniform) {
  sigmaX /= sqrt(6);
  sigmaPx /= sqrt(6);
}
if(!YUniform) {
  sigmaY /= sqrt(6);
  sigmaPy /= sqrt(6);
}
  cout << "sigmaX: "<<sigmaX<<endl;
  cout << "sigmaXp: "<<sigmaPx/pMean<<endl;
  cout << "sigmaPx: "<<sigmaPx<<endl;
  cout << "sigmaY: "<<sigmaY<<endl;
  cout << "sigmaYp: "<<sigmaPy/pMean<<endl;
  cout << "sigmaPy: "<<sigmaPy<<endl;




string outputFileName = "events.root";


TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
outputFile->mkdir("VirtualDetector");
outputFile->cd("VirtualDetector");
TNtuple *outputNtuple = new TNtuple("tracklist","Track List","x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight");
  
Float_t z=zStart;
Float_t Pz=pMean;

TRandom2 r(0);

for(int i=0;i<nEvents;i++) {
// Generate uniform distribution in X
  Float_t x,Px;  
  if(XUniform) {
   do {
     x = 2*(.5-r.Uniform());
     Px = 2*(.5-r.Uniform());
   } while(quadAdd(x,Px)>1.);  
  } else {
   r.Rannor(x,Px);
  }
  
  x *= sigmaX;
  Px *= sigmaPx;
  
  Float_t y; 
  Float_t Py; 

  if(YUniform) {
   do {
     y = 2*(.5-r.Uniform());
     Py = 2*(.5-r.Uniform());
   } while(quadAdd(y,Py)>1.);  
  } else {
   r.Rannor(y,Py);
  }
   
  y *= sigmaY;
  Py *= sigmaPy;
  
  Float_t PDGid = 2212.;
  Float_t EventID = (Float_t) i+1;
  Float_t ParentID = 0.;
  Float_t Weight = 1.;
  outputNtuple->Fill(x,y,z,Px,Py,Pz,0.,PDGid,EventID,0.,ParentID,Weight); 
  
}
outputFile->Write();
outputFile->Close();

}
Float_t quadAdd(Float_t x, Float_t y) {
   return(sqrt(x*x+y*y));
}
