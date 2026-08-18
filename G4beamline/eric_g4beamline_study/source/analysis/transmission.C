//
// Program to calculate transmission efficiency
// for a file from the Mu2e beamline. It used to be much more complicated
// Now it just counts the number of particles that hit the end.
//
#include <TNtuple.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <iostream>


double targetRadius=5;  // Target radius in mm (loose cuts)

bool targetCut(double,double);


void transmission(string filename,string startNtuple="Z2201",string endNtuple="Z245000") {



// Load the track ntuple
  TFile f(filename.c_str());
  TNtuple *beamStart = (TNtuple *) f.FindObjectAny(startNtuple.c_str());
  TNtuple *beamEnd = (TNtuple *) f.FindObjectAny(endNtuple.c_str());
    

// Now make some histograms to count the particles that reach the target
  TH2F target("target","Y vs X, all particles",200,-50.,50.,200,-50.,50);      
  TH2F hit("hit","Y vs X, hit target",100,-targetRadius,targetRadius,100,-targetRadius,targetRadius);

  target.SetMarkerStyle(21);
  hit.SetMarkerStyle(21);

// Get some totals
//  cout << "Filling histogram 1..."<<endl;
//  beamStart->Draw("Px:x>>totx","PDGid==2212&&Pz>8000.");
  cout << "Filling histogram 2..."<<endl;
  beamEnd->Draw("y:x>>target","Pz>100.");
  cout << "Filling histogram 3..."<<endl;
  beamEnd->Draw("y:x>>hit","Pz>100.&&targetCut(x,y)");
  
  TCanvas *c = new TCanvas(filename.c_str());
  
  c->Divide(1,2);
  
//  c->cd(1);
//  totx.Draw();
  c->cd(1);
  target.Draw();
  c->cd(2);
  hit.Draw();  
  string printFile = filename+".png";
  
  c->Print(printFile.c_str());
  
// Extract some totals
   int numTot=0;
   if(beamStart) numTot = beamStart->GetEntries();
   int numEnd = target.GetEntries();
   int numHit = hit.GetEntries();
cout << "numTot, numEnd, numHit"<<endl;
cout <<numTot<<
", "<<numEnd<<
", "<<numHit<<endl;
      
}
// See if it hits the target
bool targetCut(double x, double y) {
    double r=sqrt(x*x+y*y);
    return(r<targetRadius);
}
