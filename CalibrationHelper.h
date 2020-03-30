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

namespace stoppingcosmicmuonselection {

  class CalibrationHelper {

  public:
    CalibrationHelper();
    CalibrationHelper(const recob::PFParticle &thisParticle, art::Event const &evt, const int &plane);
    ~CalibrationHelper();

  private:
    TH1D *h_x;
    TH2D *h_yz_neg;
    TH2D *h_yz_pos;

  };
}

#endif
