/***
  Class containing useful functions for calorimetry.

*/
#ifndef CALORIMETRY_HELPER_H
#define CALORIMETRY_HELPER_H

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
#include "TH2D.h"
#include "TMath.h"

#include "DataTypes.h"
#include "TruedEdxHelper.h"

namespace stoppingcosmicmuonselection {

  class CalorimetryHelper {

  public:
    CalorimetryHelper();
    CalorimetryHelper(const recob::PFParticle &thisParticle, art::Event const &evt, const int &plane);
    ~CalorimetryHelper();

    // Get the calorimetry from the PFParticle
    void Set(const recob::PFParticle &thisParticle, art::Event const &evt, const int &plane);

    // Check if the calorimetry is valid
    bool IsValid();

    // Order the residual range with respect to the track direction
    void OrderResRange();

    // Get hit numb
    int GetHitNumb();
    // Get dQdx
    const std::vector<double> GetdQdx();
    // Get dEdx
    const std::vector<double> GetdEdx();
    // Get Residual range
    const std::vector<double> GetResRange();
    // Get Ordered residual range
    const std::vector<double> GetResRangeOrdered();
    // Get HitX
    const std::vector<double> GetHitX();
    // Get HitY
    const std::vector<double> GetHitY();
    // Get HitZ
    const std::vector<double> GetHitZ();
    // Get HitPeakTime
    const std::vector<double> GetHitPeakTime();
    // Get lifetime correction factors
    const std::vector<double> GetCorrFactor();
    // Get drift times
    const std::vector<double> GetDriftTime();
    // Get track pitches
    const std::vector<double> GetTrackPitch();
    // Get track indeces
    const std::vector<size_t> GetHitIndex();

    // Get the lifetime correction
    double LifeTimeCorr(double &ticks, const double &T0);

    // Fill 2D histo of dQdx vs residual range for hits in a given plane
    void FillHisto_dQdxVsRR(TH2D *h_dQdxVsRR);
    void FillHisto_dQdxVsRR(TH2D *h_dQdxVsRR, const double &tp_min, const double &tp_max);

    // Fill 2D histo of dQdx vs residual range for hits in a given plane. Correct by MC lifetime
    void FillHisto_dQdxVsRR_LTCorr(TH2D *h_dQdxVsRR);
    void FillHisto_dQdxVsRR_LTCorr(TH2D *h_dQdxVsRR, const double &tp_min, const double &tp_max);

    // Fill 2D histo for dQdx/dEdx with lifetime correction. dEdx taken from MC.
    void FillHisto_dQdEVsRR_LTCorr_MC(TH2D *h_dQdEVsRR, const double &tp_min, const double &tp_max);

    // Fill 2D histo for dQdx/dEdx with lifetime correction. dEdx taken from LandauVav.
    void FillHisto_dQdEVsRR_LTCorr_LV(TH2D *h_dQdEVsRR, const double &tp_min, const double &tp_max);

    // Set the parameters from the FHICL file
    void reconfigure(fhicl::ParameterSet const &p);

    // Return factor to correct dQdx values from bug in calorimetry module.
    double DqdxCorrection();

    // Reset
    void Reset();


  private:
    std::vector<anab::Calorimetry> _calos;
    int _trackHitNumb;
    std::vector<double> _dqdx;
    std::vector<double> _dedx;
    std::vector<double> _resrange;
    std::vector<double> _resrange_ord;
    std::vector<double> _hitx;
    std::vector<double> _hity;
    std::vector<double> _hitz;
    std::vector<double> _hitPeakTime;
    std::vector<double> _corr_factors;
    std::vector<double> _drift_time;
    std::vector<double> _track_pitch;
    std::vector<size_t> _hitIndex;

    bool _isValid = false;
    bool _isCalorimetrySet = false;
    bool _isData = false;

    std::string fTrackerTag;
    std::string fCalorimetryTag;
    std::string fPFParticleTag;

    TruedEdxHelper truedEdxHelper;

    protoana::ProtoDUNETrackUtils        trackUtil;
    protoana::ProtoDUNEPFParticleUtils   pfpUtil;

    const geo::GeometryCore *geom = lar::providerFrom<geo::Geometry>();
    const detinfo::DetectorProperties *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  };
}

#endif
