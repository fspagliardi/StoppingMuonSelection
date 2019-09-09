/***
  Class containing useful functions for geometry.

*/
#ifndef TRUEDEDX_HELPER_H
#define TRUEDEDX_HELPER_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "TVector3.h"
#include "TMath.h"

namespace stoppingcosmicmuonselection {

  class TruedEdxHelper {

  public:
    TruedEdxHelper();
    ~TruedEdxHelper();

    // Return MPV of dEdx according to landau-vavilov
    double LandauVav(double *x, double *p);

    // Return MPV of dEdx according to landau-vavilov
    double LandauVav(double &resRange);

    // Get relativist beta given the kinetic energy
    double GetBetaFromkEnergy(const double &kEnergy);

    // Get beta-gamma relativist factor given the kinetic energy
    double GetBetaGammaFromKEnergy(const double &kEnergy);

    // Get kinetic energy from res range. (Wrap for the spline)
    double ResRangeToKEnergy(const double &resRange);

    // Definition of the spline to go from res range to kinetic energy for a muon
    double Spline3(const double &x);


  private:
    double m_muon = 105.6583745; //MeV
    double c = 299792458;
    double LArdensity = 1.3954; // g/cm3
  };
}

#endif
