/***
  Class containing algorithms for stopping muons selection.

*/
#ifndef STOPPING_MUON_SELECTION_ALG_H
#define STOPPING_MUON_SELECTION_ALG_H

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "larsim/Simulation/LArVoxelData.h"
#include "larsim/Simulation/LArVoxelList.h"
#include "larsim/Simulation/SimListUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "GeometryHelper.h"
#include "SpacePointAlg.h"
#include "DataTypes.h"

namespace stoppingcosmicmuonselection {

  class StoppingMuonSelectionAlg {

  public:
    StoppingMuonSelectionAlg();
    ~StoppingMuonSelectionAlg();

    // See if there is track associated to this PFParticle
    bool IsPFParticleATrack(art::Event const &evt,recob::PFParticle const &thisParticle);

    // Read parameters from FHICL file
    void reconfigure(fhicl::ParameterSet const &p);

    // Determine if the PFParticle is a selected cathode crosser
    bool IsStoppingCathodeCrosser(art::Event const &evt, recob::PFParticle const &thisParticle);

    // For MC events, check if the track is associated to a cosmic track
    bool IsTrackMatchedToTrueCosmicTrack(art::Event const &evt, recob::PFParticle const &thisParticle);

    // Order reco start and end point based on Y position
    void OrderRecoStartEnd(TVector3 &start, TVector3 &end);

    // Get track from PFParticle
    const recob::Track GetTrackFromPFParticle(art::Event const &evt, recob::PFParticle const &thisParticle);

    // Determine min and max hit peak time for this track
    void SetMinAndMaxHitPeakTime(art::Event const &evt, recob::PFParticle const &thisParticle, double &minHitPeakTime, double &maxHitPeakTime);

  private:
    double INV_DBL = -9999999;
    bool _isACathodeCrosser = false;

    // Reconstructed information
    size_t _evNumber;
    double _trackT0;
    TVector3 _recoStartPoint, _recoEndPoint;
    double _theta_xz, _theta_yz;
    double _minHitPeakTime, _maxHitPeakTime;
    double _trackLength;
    double _trackID;

    // Truth information
    int _pdg;
    TVector3 _trueStartPoint, _trueEndPoint;
    double _trueStartT, _trueEndT;
    double _trueTrackID;

    // Helpers and algorithms
    GeometryHelper geoHelper;
    SpacePointAlg  spAlg;

    // Declare handle for particle inventory service
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Declare analysis utils
    protoana::ProtoDUNETruthUtils        truthUtil;
    protoana::ProtoDUNETrackUtils        trackUtil;
    protoana::ProtoDUNEPFParticleUtils   pfpUtil;

    // Parameters from FHICL
    std::string fTrackerTag;
    std::string fPFParticleTag;
          // Define cuts
    double length_cutoff_CC;
    double offsetFiducialBounds_CC;
    double thicknessStartVolume_CC;
    double cutMinHitPeakTime_CC;
    double cutMaxHitPeakTime_CC;
    double radiusBrokenTracksSearch_CC;
    double cutCosAngleBrokenTracks_CC;
    double cutCosAngleAlignment_CC;
    double cutContourAPA_CC;

  };
}

#endif
