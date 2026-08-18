//
//
// Generate tracks based on a normalized Guassian emittance, beta, and alpha in each plane.
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
void ellipse(Float_t &, Float_t &, double, double, bool uniform=false, double aCut=0.);

TRandom2 r(0);


// If emittance is negative, then it generates a uniform emittance.  Otherwise, it
// generates a Gaussian emittance

void EmitGen(string outputFileName, int nEvents, double emitX, double betaX, 
   double alphaX, double emitY, double betaY, double alphaY, double zStart=0.,
   double minEmitX=0., double minEmitY=0. ) {
// determine the operating mode
bool uniformX = (emitX<0);
bool uniformY = (emitY<0);
emitX = fabs(emitX);
emitY = fabs(emitY);


double mp=938.2726218;
double K=8000;
double E=K+mp;
double pMean = sqrt(E*E-mp*mp);
// Overwrite with Francois' number
double betagamma = pMean/mp;
  double sigmaX=sqrt(betaX*emitX/betagamma);
  double gammaX = (1+alphaX*alphaX)/betaX;
  double sigmaPx = pMean*sqrt(gammaX*emitX/betagamma)/1000.;
  double sigmaY = sqrt(betaY*emitY/betagamma);
  double gammaY = (1+alphaY*alphaY)/betaY;
  double sigmaPy = pMean*sqrt(gammaY*emitY/betagamma)/1000.;
  
cout << "E: "<<E<<" MeV"<<endl;
cout << "pMean:"<<pMean<<" MeV/c"<<endl;
cout << "SigmaX: "<<sigmaX<<" mm"<<endl;
cout << "SigmaY: "<<sigmaY<<" mm"<<endl;
cout << "SigmaPx: "<<sigmaPx<<" MeV/c"<<endl;
cout << "SigmaPy: "<<sigmaPy<<" MeV/c"<<endl;



//calculate gamma functions

Float_t x,xP,y,yP,Px,Py;







TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
outputFile->mkdir("VirtualDetector");
outputFile->cd("VirtualDetector");
TNtuple *outputNtuple = new TNtuple("tracklist","Track List","x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight");
  
Float_t z=zStart;
Float_t Pz=pMean;


for(int i=0;i<nEvents;i++) {


  double aCutX=minEmitX/emitX;
  double aCutY=minEmitY/emitY;
  ellipse(x,xP,betaX,alphaX,uniformX,aCutX);
  ellipse(y,yP,betaY,alphaY,uniformY,aCutY);
  
  x*=sqrt(emitX/betagamma);
  xP *= sqrt(emitX/betagamma)/1000.;
  
  y*=sqrt(emitY/betagamma);
  yP *= sqrt(emitY/betagamma)/1000.;
  

  Px=xP*pMean;
  Py=yP*pMean;

  
  Float_t PDGid = 2212.;
  Float_t EventID = (Float_t) i+1;
  Float_t ParentID = 0.;
  Float_t Weight = 1.;
  
   
  outputNtuple->Fill(x,y,z,Px,Py,Pz,0.,PDGid,EventID,0.,ParentID,Weight); 
  
}
outputFile->Write();
outputFile->Close();

}
//
// This generates an ellipse with unit emittance, based on beta and alpha, either
//  uniform and Gaussian.  Optionally, it can also have a minimum amplitude to study
//  large amplitude tails.  aCut = emit_min/emit
//
void ellipse(Float_t& x,Float_t& xP, double beta, double alpha, bool uniform, double aCut) {

   double amp;
   
// If in "uniform" mode, then uniformly populate a circle of unit radius
   if (uniform) {
     if(aCut>=1.) {
       cerr << "Error: Uniform emittance cut greater than emittance"<<endl;
       return;
     }
     do {
       x  = 2*(.5-r.Uniform());
       xP = 2*(.5-r.Uniform());
       amp = x*x+xP*xP;
     }while((amp>1.)||(amp<aCut));
     amp = sqrt(amp);

// Else generate a random Guassian with an emittance of 1.
   } else {

     // 2D gaussian with a minimum emittance cut
     amp = -log(r.Uniform(0.,exp(-aCut/2)));
     amp = sqrt(2*amp);
   }      

   // Distribute according to the amplitude and a random phase angle
     double ang = r.Uniform()*2*3.1415927;
     x  =  amp*sqrt(beta)*cos(ang);
     xP =  amp*sin(ang)/sqrt(beta) - alpha*x/beta;
}

Float_t quadAdd(Float_t x, Float_t y) {
   return(sqrt(x*x+y*y));
}
