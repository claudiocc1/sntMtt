#include <iostream>
#include <math.h>
#include "TMatrixD.h"
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"
#include <utility>      // for std::pair
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;


//---------------------------------------------------------------------
// This is "stolen" from the standard FastMTT code
//---------------------------------------------------------------------
std::pair<double, double> metTF(LorentzVector metP4,
                         LorentzVector nuP4,
                         TMatrixD covMET) {

  const double  aMETx = metP4.X();
  const double  aMETy = metP4.Y();

  double invCovMETxx = covMET(1,1);
  double invCovMETxy = -covMET(0,1);
  double invCovMETyx = -covMET(1,0);
  double invCovMETyy = covMET(0,0);
  double covDet = invCovMETxx*invCovMETyy - invCovMETxy*invCovMETyx;

  if( std::abs(covDet)<1E-10){
    std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!"
              <<"METx: "<<aMETx<<" METy: "<<aMETy
              << "MET covXX:" << covMET(0,0)<<" MET covXY:"<<covMET(1,0)<<" MET covYY:"<<covMET(1,1)<<std::endl;
    return std::make_pair(0,0);
  }
  double const_MET = 1./(2.*M_PI*TMath::Sqrt(covDet));
  double residualX = aMETx - (nuP4.X());
  double residualY = aMETy - (nuP4.Y());

  double pull2 = residualX*(invCovMETxx*residualX + invCovMETxy*residualY) +
    residualY*(invCovMETyx*residualX + invCovMETyy*residualY);
  pull2/=covDet;

  return std::make_pair(const_MET*TMath::Exp(-0.5*pull2), pull2);
}

//----------------------------------------------------------
// Returns x1 and x2 such that the two taus vectors are
// tau1P4 = vis1P4/x1
// tau1P4 = vis2P4/x2
//
// The inputs are the
// - lorentz vectors for met, diphotons, and the 1 visible taus
// - the elements of the MET covariance matrix
// - two booleans for the decay modes (true if hadronic)
// - the number of jets
// - an algorithm flag
//   = 1 if "like FastMtt, but with matrix element"
//   = 2 if "SnT likelihood, without Pt constraint...like Atlas"
//   = 3 if "SNT likelihood, with    Pt constraint...like Atlas"
//----------------------------------------------------------
std::pair<double, double>sntMtt(
	LorentzVector metP4,
	LorentzVector diphoP4,
	LorentzVector vis1P4,
	LorentzVector vis2P4,
	double MET_covxx,
	double MET_covyy,
	double MET_covxy,
	bool decay1,
	bool decay2,
	int njet,
	int algorithm) {

  // Check that the algorithm is valid
  if (algorithm<1 || algorithm > 3) {
    std::cerr << "Invalid algorithm "<< algorithm << std::endl;
    return std::make_pair(0., 0.);
  }

  // The likelihood is multiplied by a factor of (1/m)**(-power).
  // This is equivalent to (x1*x2)**(power/2).
  // The variable b is power/2 (in principle algorithm dependent)
  double b[3] = {1.5, 1.5, 1.5};  

  // The MET covariance matrix
  TMatrixD covMET(2, 2);
  covMET[0][0] = MET_covxx;
  covMET[0][1] = MET_covxy;
  covMET[1][0] = MET_covxy;
  covMET[1][1] = MET_covyy;
  
  // The coefficients of the (unnormalized) pdfs for x.
  double lepCoeff[5];
  double hadCoeff[5];
  lepCoeff[0] =   6.4627;
  lepCoeff[1] =   7.11987;
  lepCoeff[2] = -24.7644;
  lepCoeff[3] =  11.1894;
  lepCoeff[4] =   0;
  hadCoeff[0] =  -20.8584; 
  hadCoeff[1] =  -103.46;
  hadCoeff[2] =  28912.2;
  hadCoeff[3] =  -40754;
  hadCoeff[4] =  16937.9;
  double *coeff1 = lepCoeff;
  double *coeff2 = lepCoeff;
  if (decay1) coeff1 = hadCoeff;
  if (decay2) coeff2 = hadCoeff;

  // The parameters of the (unnormalized) PT pdf
  double p0, p1, p2;
  if (njet == 0) {
    p0 = 4.68036e+01;
    p1 = 4.07267e+02; 
    p2 = 2.52178e+01;
  } else {
    // p0 = 1.27005e+02;
    // p1 = 5.61223e+02;
    // p2 = 1.58225e+01;
    p0 = 4.02985e+01;
    p1 = 3.11476e+02;
    p2 = 6.77842;
  }
  
  
  // scan over x1 and x2
  double lik, x1, x2, bestLik, bestx1, bestx2;
  bestLik = -9999;
  bestx1 = -1;
  bestx2 = -1;
  int nGridPoints     = 100;
  double gridFactor   = 1./nGridPoints;    
  LorentzVector visP4     = vis1P4 + vis2P4;
    
  for(int ix1 = 1; ix1<nGridPoints; ++ix1){
    x1      = ix1*gridFactor;
    LorentzVector test1P4   = vis1P4/x1;
    
    for(int ix2 = 1; ix2<nGridPoints;++ix2){
      x2      = ix2*gridFactor;
      LorentzVector test2P4   = vis2P4/x2;
      LorentzVector testDiTau = test1P4 + test2P4;
      LorentzVector testMET   = test1P4 + test2P4 - visP4;
      LorentzVector testP4HH  = test1P4 + test2P4 + diphoP4;
      double        pTHH      = testP4HH.Pt();


      // The FastMTT-like piece
      if (algorithm == 1) {	
	// stuff associated with the limits of integration.
	//fastMTT has an additional factor of 1.15 for the mVS2 variable
	double mVS2     = visP4.M()*visP4.M() /(testDiTau.M()*testDiTau.M());
	double alpha    = mVS2;
	double mVisLeg1 = vis1P4.M();
	double mVisLeg2 = vis2P4.M();
	double mTau     = 1.776;
	double mVis1OverTauSquare = std::pow(mVisLeg1/mTau, 2);
	double mVis2OverTauSquare = std::pow(mVisLeg2/mTau, 2);
	double x1Min = std::min(1.0, mVis1OverTauSquare);
	double x2Min = std::max(mVis2OverTauSquare, mVS2);
	double x2Max = std::min(1.0, mVS2/x1Min);
	if(x2Max<x2Min) {
	  lik = 0.0;
	  // std::cout << x2Max << " " << x2Min << std::endl;
	  continue;
	}

	// The integral
	lik = 0;
	for (int j=0; j<5; j++) {
	  double temp = 0.;
	  for (int k=0; k<5; k++) {
	    if (k == j) {
	      temp = temp + coeff2[k] * log(x2Max/x2Min);
	    } else {
	      int r = k-j;
	      temp = temp + coeff2[k] * (pow(x2Max,r) - pow(x2Min,r))/r;
	    }
	  }
	  lik = lik + coeff1[j]*pow(alpha,j)*temp;   
	}	
	lik = mVS2*lik/testDiTau.M();
      

      // The Snt-like piece 
      } else {
	double lik1=0.;
	double lik2=0.;
	for (int i=0; i<5; i++) {
	  lik1 = lik1 + coeff1[i]*pow(x1,i);
	  lik2 = lik2 + coeff2[i]*pow(x2,i);
	}
	lik = lik1*lik2;
	if (algorithm == 3) {
	  lik = lik * pTHH*pow((1+pTHH/p1), -p2);
	}
      }

      // The MET_TF piece is common to every algorithm
      lik = lik * metTF(metP4, testMET, covMET).first;

      // The fudge factor, aka, the "regularization factor".
      lik = lik * pow(x2*x1, b[algorithm-1]);
		  
      // if this is the maximum or the 1st time, store it
      if ( ix1==1 && ix2==1 ) {
        bestLik    = lik;
        bestx1     = x1;
        bestx2     = x2;
      } else if ( lik > bestLik) {
        bestLik    = lik;
        bestx1     = x1;
        bestx2     = x2;
      }
    } // close iX2 loop
  }   // close iX1 loop
  return std::make_pair(bestx1, bestx2);
}

