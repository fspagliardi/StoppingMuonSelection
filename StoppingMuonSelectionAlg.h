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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "GeometryHelper.h"
#include "HitHelper.h"
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

    // Determine if the PFParticle is a selected anode crosser
    bool IsStoppingAnodeCrosser(art::Event const &evt, recob::PFParticle const &thisParticle);

    // Work out t0 for anode crossers.
    double CorrectPosAndGetT0(TVector3 &_recoStartPoint, TVector3 &_recoEndPoint);

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

    // Set MCParticle properties
    void SetMCParticleProperties(art::Event const &evt, recob::PFParticle const &thisParticle);

    // Check if the true track associated with that PFParticle is a stopping muon
    bool IsTrueParticleAnAnodeCrossingStoppingMuon(art::Event const &evt, recob::PFParticle const &thisParticle);

    // Check if the true track associated with that PFParticle is a stopping muon
    bool IsTrueParticleACathodeCrossingStoppingMuon(art::Event const &evt, recob::PFParticle const &thisParticle);

    // N-1 cuts for Cathode crossers
    bool NMinus1Cathode(const std::string &excludeCut, art::Event const &evt, const recob::PFParticle &thisParticle);

    // N-1 cuts for Anode crossers
    bool NMinus1Anode(const std::string &excludeCut, art::Event const &evt, const recob::PFParticle &thisParticle);

    // N-1 cuts for Cathode crossers (simple version)
    bool NMinus1CathodeSimple(const std::string &excludeCut, art::Event const &evt, const recob::PFParticle &thisParticle);

    // Get the property for this track.
    const trackProperties GetTrackProperties();

    // Reset function
    void Reset();

  private:
    bool _isACathodeCrosser = false;
    bool _isAnAnodeCrosser = false;
    bool _isPFParticleATrack = false;
    bool _areMCParticlePropertiesSet = false;

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

    trackProperties trackInfo;

    // Helpers and algorithms
    GeometryHelper geoHelper;
    HitHelper      hitHelper;

    // Declare handle for particle inventory service
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Declare handle for detector properties
    const detinfo::DetectorProperties *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

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
    double length_cutoff_AC;
    double offsetYStartPoint_AC;
    double offsetZStartPoint_AC;
    double cutMinHitPeakTime_AC;
    double cutMaxHitPeakTime_AC;
    double radiusBrokenTracksSearch_AC;
    double cutCosAngleBrokenTracks_AC;
    double cutCosAngleAlignment_AC;
    double cutContourAPA_AC;
    double offsetFiducialBounds_AC;

  };
}

#endif
