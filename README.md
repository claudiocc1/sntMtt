
Fitting tau tau mass in the context of HH->gamma gamma tau tau analysis

Usage:
```
    std::pair<double, double> p4_scalings = sntMtt(metP4, diphoP4, vis1P4, vis2P4, MET_covxx,
                                                   MET_covyy, MET_covxy, dec1, dec2, njet,algo);
    LorentzVector tau1_p4  = vis1_p4/p4_scalings.first;
    LorentzVector tau2_p4  = vis2_p4/p4_scalings.second;
    LorentzVector ditau_p4 = tau1_p4 + tau2_p4;
    LorentzVector HH_p4    = dipho_p4 + ditau_p4;
    double mTauTau         = ditau_p4.M()
```
The sntMTT function:
```
//----------------------------------------------------------
// Returns x1 and x2 such that the two taus vectors are
// tau1P4 = vis1P4/x1
// tau1P4 = vis2P4/x2
//
// The inputs are the
// - lorentz vectors for met, diphotons, and the visible taus
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
        int algorithm)
```
The three algorithms and their performance on signal Monte Carlo compared to FastMtt are described in the pdf file. 

Algorithms 1 and 2 have resolution similar to FastMtt.

Algorithm 3 improves the resolution by about 10-15% in the zero jet bin.


