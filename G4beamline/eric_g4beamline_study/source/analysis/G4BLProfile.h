#ifndef __G4BLPROFILE_H__
#define __G4BLPROFILE_H__
//
//  This is a quick and dirty program to convert a G4 profile to a ROOT NTUPLE
//
// Usage:
//   G4BLProfile f(string filename);
//   TNtuple *n=f.getNtuple();
// The Ntuple will be called "G4BLProfile"
//
//  2014-02-12  E. Prebys Original
//  2014-02-13  E. Prebys Changed name
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <TNtuple.h>

using namespace std;

class G4BLProfile{
  public:
    TNtuple *n;
    G4BLProfile(string); 
    TNtuple *getNtuple();
};

#endif
      