#include "THcPShHit.h"
#include "TMath.h"

#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>

using namespace std;

// Track class for the SHMS calorimeter calibration.
// Comprises the spectrometer track parameters and calorimeter hits.
//

// Container (collection) of hits and its iterator.
//
typedef vector<THcPShHit*> THcPShHitList;
typedef THcPShHitList::iterator THcPShHitIt;

class THcPShTrack {

  Double_t P;   // track momentum
  Double_t X;   // at the Preshower face
  Double_t Xp;  // slope
  Double_t Y;   // at the Preshower face
  Double_t Yp;  // slope

  THcPShHitList Hits;

 public:
  THcPShTrack();
  THcPShTrack(Double_t p, Double_t x, Double_t xp, Double_t y, Double_t yp);
  ~THcPShTrack();

  void Reset(Double_t p, Double_t x, Double_t xp, Double_t y, Double_t yp);

  void AddHit(Double_t adc, Double_t edep, UInt_t blk_number);

  THcPShHit* GetHit(UInt_t k);

  UInt_t GetNhits() {return Hits.size();};

  void Print(ostream & ostrm);

  void SetEs(Double_t* alpha);

  Double_t Enorm();

  Double_t GetP() {return P*1000.;}     //MeV

  Float_t Ycor(Double_t, Int_t);        // coord. corection for Preshower module

  // Coordinate correction constants for Preshower blocks
  //
  static const Double_t fAcor = 106.73;
  static const Double_t fBcor = 2.329;

  // Calorimeter geometry constants.
  //
  //  static const Double_t fZbl = 10;   //cm, Preshower block transverse size
  static const UInt_t fNrows_pr = 14;    //Row number for Preshower
  static const UInt_t fNrows_sh = 16;    //Row number for Shower
  static const UInt_t fNcols_pr =  2;    //2 columns in Preshower
  static const UInt_t fNcols_sh = 14;    //14 columnsin Shower
  static const UInt_t fNpmts_pr = fNrows_pr*fNcols_pr;
  static const UInt_t fNpmts = fNpmts_pr + fNrows_sh*fNcols_sh;;

};

//------------------------------------------------------------------------------

THcPShTrack::THcPShTrack() { };

THcPShTrack::THcPShTrack(Double_t p,
		       Double_t x, Double_t xp, Double_t y, Double_t yp) {
  P = p;
  X = x;
  Xp = xp;
  Y = y;
  Yp =yp;
};

//------------------------------------------------------------------------------

void THcPShTrack::Reset(Double_t p,
		       Double_t x, Double_t xp, Double_t y, Double_t yp) {

  // Reset track parameters, clear hit list.

  P = p;
  X = x;
  Xp = xp;
  Y = y;
  Yp =yp;
  Hits.clear();
};

//------------------------------------------------------------------------------

void THcPShTrack::AddHit(Double_t adc, Double_t edep, UInt_t blk_number) {

  // Add a hit to the hit list.

  THcPShHit* hit = new THcPShHit(adc, blk_number);
  hit->SetEdep(edep);
  Hits.push_back(hit);
};

//------------------------------------------------------------------------------

THcPShHit* THcPShTrack::GetHit(UInt_t k) {
  THcPShHitIt it = Hits.begin();
  for (UInt_t i=0; i<k; i++) it++;
  return *it;
}

void THcPShTrack::Print(ostream & ostrm) {

  // Output the track parameters and hit list through the stream ostrm.

  ostrm << P << " " << X << " " << Xp << " " << Y << " " << Yp << " "
	<< Hits.size() << endl;

  for (THcPShHitIt iter = Hits.begin(); iter != Hits.end(); iter++) {
    (*iter)->Print(ostrm);
  };

};

//------------------------------------------------------------------------------

THcPShTrack::~THcPShTrack() {
  for (THcPShHitIt i = Hits.begin(); i != Hits.end(); ++i) {
    delete *i;
    *i = 0;
  }
};

//------------------------------------------------------------------------------

void THcPShTrack::SetEs(Double_t* alpha) {

  // Set hit energy depositions by use of calibration (gain) constants alpha.
  
  for (THcPShHitIt iter = Hits.begin(); iter != Hits.end(); iter++) {
  
    Double_t adc = (*iter)->GetADC();
    UInt_t nblk = (*iter)->GetBlkNumber();

    if(nblk <= fNrows_pr*fNcols_pr) {  //Preshower, correct for Y coordinate
      Int_t ncol = 1;
      if (nblk > fNrows_pr) ncol = 2;
      (*iter)->SetEdep(adc*Ycor(Y,ncol)*alpha[nblk-1]);
    }
    else                               //Shower
      (*iter)->SetEdep(adc*alpha[nblk-1]);

  };

}

//------------------------------------------------------------------------------

Double_t THcPShTrack::Enorm() {

  // Normalized to the track momentum energy depostion in the calorimeter.

  Double_t sum = 0;

  for (THcPShHitIt iter = Hits.begin(); iter != Hits.end(); iter++) {
    sum += (*iter)->GetEdep();
  };

  return sum/P/1000.;
}

//------------------------------------------------------------------------------

// Coordinate correction for Preshower modules.
// Fit to GEANT pion data @ 5 GeV/c (Simon).

Float_t THcPShTrack::Ycor(Double_t yhit, Int_t ncol) {

  // Warn if hit does not belong to Preshower.
  if (ncol > fNcols_pr || ncol < 1)
    cout << '*** THcPShTrack::Ycor: wrong ncol = ' << ncol << ' ***' << endl;

  // Check hit coordinate with fired block column.
  if ((yhit < 0. && ncol == 1) || (yhit > 0. && ncol == 2))
    return 1./(1. + TMath::Power(TMath::Abs(yhit)/fAcor, fBcor));
  else
    return 1.;

}
