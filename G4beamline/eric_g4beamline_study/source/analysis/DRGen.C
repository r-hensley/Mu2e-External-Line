//
//
// Generate tracks based on DR extraction model and MARS
//  File is currently hardcoded to use ExtTracks.root as input
//
//  Dec. 11, 2015 EjP  Converted from old (cm, rad) to new (mm, mrad) format
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


// Ntuple created from DR tracks
#define ExtTracks_cxx
#include "ExtTracks.h"

void ExtTracks::Loop() {};   // expects to see this


using namespace std;

TRandom2 r(0);


// If emittance is negative, then it generates a uniform emittance.  Otherwise, it
// generates a Gaussian emittance

void DRGen(string outputFileName ,double minXEmittance=0., double minYEmittance=0.,int nSmear=1, bool rotateY=true) {
// WARNING!!! These cuts are hardcoded.  This is the range of input X which
// defines the core (ie, not the scattered tail).  It may change if the 
// distribution changes.
double xMin=-.0025;
double xMax=.0025;

// These are also hardcoded to Eliana's optics!
//
// Starting point and optics, from Eliana's optics file
double zStart = 2200.3;  //Downstream end of C magnet, in mm
double betaX = 7.098112;
double alphaX = 1.10034;
double gammaX = (1+alphaX*alphaX)/betaX;
double betaY = 9.400828;
double alphaY = -1.449754;
double gammaY = (1+alphaY*alphaY)/betaY;




//Smearing factors
double sigmaSmearX=.05;
double sigmaSmearPx=.1;
double sigmaSmearY=.1;
double sigmaSmearPy=.2;

//calculate gamma functions


double mp=938.272;


TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
outputFile->mkdir("VirtualDetector");
outputFile->cd("VirtualDetector");

TNtuple *diagNtuple = new TNtuple("diagnostics","Diagnostic Ntuple","old_y:old_Py:normX:normY:theta:A:y:Py");
TNtuple *outputNtuple = new TNtuple("tracklist","Track List","x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight");
Float_t x,y,z,Px,Py,Pz,t,PDGid,EventID,TrackID,ParentID,Weight;
// These things will be the same for every particle
z = zStart;       // End of C magnet
PDGid = 2212.;    // proton
TrackID = 1.; 
ParentID = 0.;
Weight = 1.;     
t = 0.;


// Loop over the entries in the MARS ntuple

ExtTracks *ext = new ExtTracks();   // This will automatically load the tree from
                                    // ExtTracks.root

// Loop over the entries
Long64_t nentries = ext->fChain->GetEntriesFast();

// Loop over entries to get the means, then loop again to calculate the
// optical functions.
double xSum=0;
double x2Sum=0;
double xpSum=0;
double xp2Sum=0;
double ySum=0;
double y2Sum=0;
double ypSum=0;
double yp2Sum=0;
double xxpSum=0;
double yypSum=0;

int nsum=0;
// Calculate the means and convert to meters and radians
for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = ext->LoadTree(jentry);
      if (ientry < 0) break;
      ext->fChain->GetEntry(jentry);  
      double x = ext->x/1000;
      double xp = ext->xp/1000;
      double y = ext->y/1000;
      double yp = ext->yp/1000;
      if(x<xMin||x>xMax) continue;
      xSum += x;
      x2Sum += x*x;
      xpSum += xp;
      xp2Sum += xp*xp;
      xxpSum += x*xp;
      ySum += y;
      y2Sum += y*y;
      ypSum += yp;
      yp2Sum += yp*yp;
      yypSum += y*yp;
      nsum++;
}

// Means
double xMean = xSum/nsum;
double xpMean = xpSum/nsum;
double yMean = ySum/nsum;
double ypMean = ypSum/nsum;
// Sigmas
double xSigma = sqrt(x2Sum/nsum-xMean*xMean);
double ySigma = sqrt(y2Sum/nsum-yMean*yMean);    
double xpSigma = sqrt(xp2Sum/nsum-xpMean*xpMean);
double ypSigma = sqrt(yp2Sum/nsum-ypMean*ypMean);  
// Cross correlation
double xxpCorel = xxpSum/nsum-xMean*xpMean;
double yypCorel = yypSum/nsum-yMean*ypMean;

// Calculate Optical functions
double oldEmittanceY=sqrt(ySigma*ySigma*ypSigma*ypSigma-yypCorel*yypCorel);
double oldBetaY = ySigma*ySigma/oldEmittanceY;
double oldAlphaY = -yypCorel/oldEmittanceY;
double oldGammaY = ypSigma*ypSigma/oldEmittanceY;
double GammaCheck = (1+oldAlphaY*oldAlphaY)/oldBetaY;

cout << "Calculated functions..."<<endl;
cout << "Entries processed: "<<nsum<<"/"<<nentries<<endl;
cout << "xMean: "<<xMean<<" m"<<endl;
cout << "xSigma: "<<xSigma<<" m"<<endl;
cout << "xpMean: "<<xpMean<<endl;
cout << "xpSigma: "<<xpSigma<<endl;
cout << "xxpCorel: "<<xxpCorel<<" m"<<endl;
cout << "yMean: "<<yMean<<" m"<<endl;
cout << "ySigma: "<<ySigma<<" m"<<endl;
cout << "ypMean: "<<ypMean<<endl;
cout << "ypSigma: "<<ypSigma<<endl;
cout << "yypCorel: "<<yypCorel<<" m"<<endl;

cout << "Y Emittance: "<<oldEmittanceY<<" m"<<endl;
cout << "Y Beta: "<<oldBetaY<<" m"<<endl;
cout << "Y Alpha: "<<oldAlphaY<<endl;
cout << "Y Gamma: "<<oldGammaY<<" ("<<GammaCheck<<") /m"<<endl;


      
long nTotal = 0;
long nSaved = 0;
ySum = y2Sum = ypSum = yp2Sum = yypSum = nsum = 0.;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = ext->LoadTree(jentry);
      if (ientry < 0) break;
      ext->fChain->GetEntry(jentry);  

// Calculate the momentum
      double K = ext->k*1000.;  // GeV->Mev
      double E=K+mp;
      double Ptot = sqrt(E*E-mp*mp);
      double betagamma=Ptot/mp;

      x = (ext->x-xMean*1000.);  
      
      double xp = ext->xp/1000.-xpMean;
      Px = Ptot*xp;
      y = (ext->y-yMean*1000.);
      double yp = ext->yp/1000.-ypMean;
      Py = Ptot*yp;
      Pz = sqrt(Ptot*Ptot - Py*Py - Px*Px);

      if(rotateY) {
// Test mode
        if(oldBetaY<0.) {
          oldBetaY=betaY;
          oldAlphaY=alphaY;
        }

        double old_y=y;
        double old_Py=Py;
        
// Do a Twiss transformation in the Y plane
// Calculate the amplitude and angle in the existing system
        double oldGammaY = (1+oldAlphaY*oldAlphaY)/oldBetaY;
        double normY = yp*sqrt(oldBetaY) + oldAlphaY*y/1000./sqrt(oldBetaY);
        double normX = y/1000./sqrt(oldBetaY);
        double theta=0.;
        if((normX*normY)!=0.) {
          theta = atan2(normY,normX);
        }
        double A = sqrt(normX*normX+normY*normY);
        
//        if(nTotal<5.) {
//           cout << "input values: "<<y<<", "<<yp<<", "<<Py<<endl;
//           cout << "old Twiss: "<<oldBetaY<<", "<<oldAlphaY<<", "<<oldGammaY<<", "<<endl;
//           cout<< "new Twiss: "<<betaY<<", "<<alphaY<<", "<<gammaY<<", "<<endl;
//           cout <<" normX,normY,theta,A: "<<normX<<", "<<normY<<", "<<theta<<", "<<A<<endl;
//        }
// Now rotate it into the new alpha
        y = A*cos(theta);
        yp = A*sin(theta) - alphaY*y;
      
        y *= sqrt(betaY)*1000.; //convert back to mm
        yp /= sqrt(betaY);
      
        Py = Ptot*yp;
        
        diagNtuple->Fill(old_y,old_Py,normX,normY,theta,A,y,Py);
        
//        if(nTotal<5.) {
//           cout << "output values: "<<y<<", "<<yp<<", "<<Py<<endl;
//        }  
      }

// Save the initial values, so they can be smeared.      
      double x0 = x;
      double Px0 = Px;
      double y0 = y;
      double Py0 = Py;
      
      for (int is=0; is<nSmear;is++) {
      
// Calculate the emittances
        double xEmit = gammaX*x*x*1e-6 + 2*alphaX*xp*x*1e-3 + betaX*xp*xp;
        xEmit *= 1e6*betagamma;  // normalized in mm-mr
        double yEmit = gammaY*y*y*1e-6 + 2*alphaY*yp*y*1e-3 + betaY*yp*yp;
        yEmit *= 1e6*betagamma;  // normalized in mm-mr

         
// Only write it out if it's greater than the cut
        if((xEmit>minXEmittance)||(yEmit>minYEmittance)) {      
           nSaved++;
           EventID = nSaved;
           outputNtuple->Fill(x,y,z,Px,Py,Pz,0.,PDGid,EventID,TrackID,ParentID,Weight); 
        }
// Smear for the next pass
        x = x0 + sigmaSmearX*r.Gaus();
        Px = Px0 + sigmaSmearPx*r.Gaus();
        y = y0 + sigmaSmearY*r.Gaus();
        Py = Py0 + sigmaSmearPy*r.Gaus();
        
      }   
      
// Calculate the output optical functions
      if((x>=(xMin*1000.))&&(x<=(xMax*1000.))) {
        ySum += y0/1000.;
        y2Sum += y0*y0/1000000.;
        ypSum += Py0/Pz;
        yp2Sum += Py0*Py0/(Pz*Pz);
        yypSum += y0*Py0/Pz/1000.;
        nsum++;
      }      
      nTotal++;
  }

cout << "Processed "<<nTotal<<" tracks, saved "<<nSaved<<".  Closing file..."<<endl;
// Means
yMean = ySum/nsum;
ypMean = ypSum/nsum;
// Sigmas
xpSigma = sqrt(xp2Sum/nsum-xpMean*xpMean);
ypSigma = sqrt(yp2Sum/nsum-ypMean*ypMean);  
// Cross correlation
yypCorel = yypSum/nsum-yMean*ypMean;

// Calculate Optical functions
double newEmittanceY=sqrt(ySigma*ySigma*ypSigma*ypSigma-yypCorel*yypCorel);
double newBetaY = ySigma*ySigma/newEmittanceY;
double newAlphaY = -yypCorel/newEmittanceY;
double newGammaY = ypSigma*ypSigma/newEmittanceY;
GammaCheck = (1+newAlphaY*newAlphaY)/newBetaY;

cout << "Calculated output optics..."<<endl;
cout << "nsum: "<<nsum<<"/"<<nTotal<<endl;
cout << "Y Emittance: "<<newEmittanceY<<" ("<<oldEmittanceY<<") m"<<endl;
cout << "Y Beta: "<<newBetaY<<" ("<<betaY<<") m"<<endl;
cout << "Y Alpha: "<<newAlphaY<<" ("<<alphaY<<")"<<endl;
cout << "Y Gamma: "<<newGammaY<<" ("<<GammaCheck<<") /m"<<endl;


outputFile->Write();
outputFile->Close();

cout << "Done."<<endl;

}


