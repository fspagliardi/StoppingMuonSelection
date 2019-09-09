/***
  Class containing useful functions for geometry.

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

namespace stoppingcosmicmuonselection {

  class CalorimetryHelper {

  public:
    CalorimetryHelper();
    CalorimetryHelper(const recob::PFParticle &thisParticle, art::Event const &evt);
    ~CalorimetryHelper();

    // Get the calorimetry from the PFParticle
    void Set(const recob::PFParticle &thisParticle, art::Event const &evt);

    // Check if the calorimetry is valid
    bool IsValid();

    // Order the residual range with respect to the track direction
    void OrderResRange();

    // Get hit numb
    int GetHitNumb(const int &planeNumb);
    // Get dQdx
    double *GetdQdx(const int &planeNumb);
    // Get dEdx
    double *GetdEdx(const int &planeNumb);
    // Get Residual range
    double *GetResRange(const int &planeNumb);
    // Get Ordered residual range
    double *GetResRangeOrdered(const int &planeNumb);
    // Get HitX
    double *GetHitX(const int &planeNumb);
    // Get HitY
    double *GetHitY(const int &planeNumb);
    // Get HitZ
    double *GetHitZ(const int &planeNumb);
    // Get HitPeakTime
    double *GetHitPeakTime(const int &planeNumb);
    // Get lifetime correction factors
    double *GetCorrFactor(const int &planeNumb);
    // Get drift times
    double *GetDriftTime(const int &planeNumb);
    // Get track pitches
    double *GetTrackPitch(const int &planeNumb);

    // Get the lifetime correction
    double LifeTimeCorr(double &ticks, const double &T0);

    // Get 2D histo of dQdx vs residual range for hits in a given plane
    TH2D *GetHisto_dQdxVsRR(const int &planeNumb);

    TH2D *GetHisto_dQdxVsRR(const int &planeNumb, const double &tp_min, const double &tp_max);


    // Set the parameters from the FHICL file
    void reconfigure(fhicl::ParameterSet const &p);

    // Reset
    void Reset();


  private:
    std::vector<anab::Calorimetry> _calos;
    int _trackHitNumb[3];
    double _dqdx[3][3000];
    double _dedx[3][3000];
    double _resrange[3][3000];
    double _resrange_ord[3][3000];
    double _hitx[3][3000];
    double _hity[3][3000];
    double _hitz[3][3000];
    double _hitPeakTime[3][3000];
    double _corr_factors[3][3000];
    double _drift_time[3][3000];
    double _track_pitch[3][3000];

    bool _isValid = false;
    bool _isCalorimetrySet = false;
    bool _isData = false;

    std::string fTrackerTag;
    std::string fCalorimetryTag;
    std::string fPFParticleTag;

    protoana::ProtoDUNETrackUtils        trackUtil;
    protoana::ProtoDUNEPFParticleUtils   pfpUtil;

    const geo::GeometryCore *geom = lar::providerFrom<geo::Geometry>();
    const detinfo::DetectorProperties *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  };
}

#endif
