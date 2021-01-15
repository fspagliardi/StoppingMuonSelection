/***
  Class containing useful functions for geometry.

*/
#ifndef TRUEDEDX_HELPER_CXX
#define TRUEDEDX_HELPER_CXX

#include "TruedEdxHelper.h"

namespace stoppingcosmicmuonselection {

  TruedEdxHelper::TruedEdxHelper() {
    TFile inputFile("./MCdEdxSuperBinning.root");
    _h_dEdx = (TH1D *)inputFile.Get("h_MPV");
  }

  TruedEdxHelper::~TruedEdxHelper() {

  }

  // Return MPV of dEdx according to landau-vavilov
  double TruedEdxHelper::LandauVav(double *x, double *p, const double &LArdensity)  {
    double resRange = x[0];
    double X = p[0] * LArdensity;
  	double d = p[1];
  	double Z = 18;
    double mc2 = 0.5109989461; // MeV
    double K = 0.307075; // MeV cm2 / mol
    double A = 39.948; // g/mol
    double z = 1;
    double I = 0.000188;
    double j = 0.200;
    double kEn = ResRangeToKEnergy(resRange);
  	//double en = x[0];
    double yb = GetBetaGammaFromKEnergy(kEn);
    double b = GetBetaFromkEnergy(kEn);
    double E = (K/2) * (Z/A) * z * z * (X / TMath::Power(b,2) );
    double MPV = LArdensity * 1/X * E *(TMath::Log(2*mc2*yb*yb/I)+TMath::Log(E/I)+j-(b*b)-d);
    return MPV;
  }

  // Return MPV of dEdx according to landau-vavilov
  double TruedEdxHelper::LandauVav(double &resRange, const double &trackPitch, const double &LArdensity) {
    double kEn = ResRangeToKEnergy(resRange);
    double yb = GetBetaGammaFromKEnergy(kEn);
    double d = DensityEffect(yb); // density correction
    double pars[2] = {trackPitch, d};
    return LandauVav(&resRange, pars, LArdensity);
  }

  // Get the dEdx from the MC simulation.
  double TruedEdxHelper::GetMCdEdx(const double &resRange) {
    double MCdEdx = _h_dEdx->GetBinContent(_h_dEdx->FindBin(resRange));
    return MCdEdx;
  }

  // Work out the density effect based on Sternheimer parametrization. (https://journals.aps.org/prb/pdf/10.1103/PhysRevB.3.3681)
  double TruedEdxHelper::DensityEffect(const double &yb) {
    double x = TMath::Log10(yb);
    // See also https://www.sciencedirect.com/science/article/pii/S0092640X01908617
    // page 204
    // values eventually taken from http://pdg.lbl.gov/2019/AtomicNuclearProperties/MUE/muE_liquid_argon.pdf
    double C = -5.2146; // Always negative
    double x1 = 3.0;
    double x0 = 0.2;
    double m = 3.0;
    double a = 0.19559;
    if (x >= x1) {
      return 2*TMath::Log(10)*x+C;
    }
    else if (x <= x0)
      return 0;
    else {
      return 2*TMath::Log(10)*x + C + a*TMath::Power(x1-x,m);
    }
  }

  // Get relativist beta given the kinetic energy
  double TruedEdxHelper::GetBetaFromkEnergy(const double &kEnergy)  {
    double beta = (TMath::Sqrt((kEnergy+m_muon)*(kEnergy+m_muon) - m_muon*m_muon))  / (kEnergy+m_muon);

    return beta;
  }

  // Get beta-gamma relativist factor given the kinetic energy
  double TruedEdxHelper::GetBetaGammaFromKEnergy(const double &kEnergy) {
    double yb = TMath::Sqrt((kEnergy+m_muon)*(kEnergy+m_muon) - (m_muon*m_muon)) / m_muon;

    return yb;
  }

  // Get kinetic energy from res range. (Wrap for the spline)
  double TruedEdxHelper::ResRangeToKEnergy(const double &resRange) {
    return Spline3(resRange);
  }

  // Definition of the spline to go from res range to kinetic energy for a muon
  double TruedEdxHelper::Spline3(const double &x) {
     const int fNp = 49, fKstep = 0;
     const double fDelta = -1, fXmin = 0.00202794, fXmax = 441.476;
     const double fX[49] = { 0.00202794, 0.00661748, 0.0118338, 0.020788, 0.031053,
                          0.0509814, 0.0742837, 0.10086, 0.130516, 0.163324,
                          0.198997, 0.237607, 0.279011, 0.370057, 0.471562,
                          0.583166, 0.70437, 0.974212, 1.27937, 1.79585,
                          2.37894, 3.48066, 4.72636, 6.09742, 7.5788,
                          9.15473, 10.8166, 12.5501, 14.3553, 18.1304,
                          22.0917, 26.2106, 30.4441, 39.2049, 48.2235,
                          62.0774, 76.1461, 99.8567, 123.567, 147.278,
                          170.845, 194.198, 217.407, 240.473, 263.395,
                          308.739, 353.438, 397.708, 441.476 };
     const double fY[49] = { 1, 1.2, 1.4, 1.7, 2,
                          2.5, 3, 3.5, 4, 4.5,
                          5, 5.5, 6, 7, 8,
                          9, 10, 12, 14, 17,
                          20, 25, 30, 35, 40,
                          45, 50, 55, 60, 70,
                          80, 90, 100, 120, 140,
                          170, 200, 250, 300, 350,
                          400, 450, 500, 550, 600,
                          700, 800, 900, 1000 };
     const double fB[49] = { 46.5968, 40.824, 36.2021, 31.1895, 27.5109,
                          23.0916, 20.0304, 17.7713, 16.007, 14.5855,
                          13.467, 12.4881, 11.6882, 10.3741, 9.38014,
                          8.5832, 7.94817, 6.94495, 6.22119, 5.44824,
                          4.8863, 4.24693, 3.81315, 3.49923, 3.26646,
                          3.08479, 2.94289, 2.82477, 2.72295, 2.58223,
                          2.47118, 2.39169, 2.33263, 2.2442, 2.19302,
                          2.14516, 2.12009, 2.106, 2.10846, 2.11271,
                          2.13175, 2.14843, 2.16084, 2.17432, 2.18876,
                          2.22258, 2.24882, 2.27045, 2.30055 };
     const double fC[49] = { -715.905, -541.909, -344.151, -215.653, -142.704,
                          -79.0586, -52.3073, -32.7012, -26.7907, -16.5347,
                          -14.8218, -10.5307, -8.7881, -5.6451, -4.14744,
                          -2.99332, -2.24607, -1.47174, -0.899992, -0.596609,
                          -0.367106, -0.213227, -0.134998, -0.0939646, -0.0631621,
                          -0.0521182, -0.0332642, -0.0348786, -0.0215249, -0.0157522,
                          -0.0122798, -0.00701949, -0.00693155, -0.00316226, -0.0025128,
                          -0.000942072, -0.000839419, 0.000245231, -0.000141504, 0.000320787,
                          0.000487138, 0.000226887, 0.000307886, 0.00027668, 0.000353125,
                          0.000392741, 0.000194312, 0.000294337, 43.7679 };
     const double fD[49] = { 12637.1, 12637.1, 4783.55, 2368.86, 1064.56,
                          382.671, 245.913, 66.4331, 104.203, 16.0058,
                          37.046, 14.0292, 11.507, 4.91824, 3.44704,
                          2.05509, 0.956519, 0.624541, 0.195804, 0.131198,
                          0.0465572, 0.020933, 0.00997609, 0.00693107, 0.00233595,
                          0.00378163, -0.000310421, 0.00246584, 0.000509722, 0.000292195,
                          0.0004257, 6.92434e-06, 0.000143416, 2.40044e-05, 3.77928e-05,
                          2.43218e-06, 1.52485e-05, -5.43688e-06, 6.49908e-06, 2.35285e-06,
                          -3.71483e-06, 1.16332e-06, -4.50973e-07, 1.11163e-06, 2.91229e-07,
                          -1.47974e-06, 7.53155e-07, 7.53155e-07, 20.5494 };
     int klow=0;
     if(x<=fXmin) klow=0;
     else if(x>=fXmax) klow=fNp-1;
     else {
       if(fKstep) {
         // Equidistant knots, use histogramming
         klow = int((x-fXmin)/fDelta);
         if (klow < fNp-1) klow = fNp-1;
       } else {
         int khig=fNp-1, khalf;
         // Non equidistant knots, binary search
         while(khig-klow>1)
           if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
           else khig=khalf;
       }
     }
     // Evaluate now
     double dx=x-fX[klow];
     return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
  }

} // end of namespace stoppingcosmicmuonselection

#endif
