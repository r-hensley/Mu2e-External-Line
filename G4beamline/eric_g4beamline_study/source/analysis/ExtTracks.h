//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 11 14:19:00 2015 by ROOT version 6.04/02
// from TTree ExtTracks/Tracks produced by MARS simulation of Delivery ring extraction
// found on file: ExtTracks.root
//////////////////////////////////////////////////////////

#ifndef ExtTracks_h
#define ExtTracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ExtTracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         x;
   Float_t         xp;
   Float_t         y;
   Float_t         yp;
   Float_t         t;
   Float_t         k;

   // List of branches
   TBranch        *b_x;   //!
   TBranch        *b_xp;   //!
   TBranch        *b_y;   //!
   TBranch        *b_yp;   //!
   TBranch        *b_t;   //!
   TBranch        *b_k;   //!

   ExtTracks(TTree *tree=0);
   virtual ~ExtTracks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ExtTracks_cxx
ExtTracks::ExtTracks(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ExtTracks.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ExtTracks.root");
      }
      f->GetObject("ExtTracks",tree);

   }
   Init(tree);
}

ExtTracks::~ExtTracks()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ExtTracks::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ExtTracks::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ExtTracks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("xp", &xp, &b_xp);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("yp", &yp, &b_yp);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("k", &k, &b_k);
   Notify();
}

Bool_t ExtTracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ExtTracks::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ExtTracks::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ExtTracks_cxx
