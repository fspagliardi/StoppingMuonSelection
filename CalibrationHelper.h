/***
  Class containing useful functions for calibration.

*/
#ifndef CALIBRATION_HELPER_H
#define CALIBRATION_HELPER_H

#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"

#include "DataTypes.h"
#include "SceHelper.h"

namespace stoppingcosmicmuonselection {

  class CalibrationHelper {

  public:
    CalibrationHelper();
    CalibrationHelper(const std::string &filetype);
    ~CalibrationHelper();

    // Get the histos.
    void Set(art::Event const &evt);

    // Get X correction factor.
    double GetXCorr(const TVector3 &hitPos);

    // Get YZ correction factor.
    double GetYZCorr(const TVector3 &hitPos);

    // Get both factors at the same time.
    double GetXYZCorr(const TVector3 &hitPos);

    // Get vector of factors for X.
    std::vector<double> GetXCorr_V(const std::vector<double> &hit_xs);

    // Get vector of factors for YZ.
    std::vector<double> GetYZCorr_V(const std::vector<double> &hit_xs, const std::vector<double> &hit_ys, const std::vector<double> &hit_zs);

    // Get vector of directions.
    std::vector<TVector3> GetHitDirVec(const std::vector<double> &hit_xs, const std::vector<double> &hit_ys, const std::vector<double> &hit_zs);

    // Get vector of Fields.
    std::vector<TVector3> GetHitPosField(const std::vector<double> &hit_x, const std::vector<double> &hit_y, const std::vector<double> &hit_z);

    // Get vector of angles phi.
    std::vector<double> PitchFieldAngle(const std::vector<double> &hit_xs, const std::vector<double> &hit_ys, const std::vector<double> &hit_zs);

  private:
    TH1D *h_x;
    TH2D *h_yz_neg;
    TH2D *h_yz_pos;

    SceHelper sceHelper;

  };
}

#endif
