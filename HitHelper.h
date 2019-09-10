/***
  Class containing useful functions for hit.

*/
#ifndef HIT_HELPER_H
#define HIT_HELPER_H

#include "lardataobj/Simulation/SimChannel.h"
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "fhiclcpp/ParameterSet.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "TVector3.h"
#include "TH2D.h"
#include "TMath.h"

#include "DataTypes.h"

namespace stoppingcosmicmuonselection {

  class HitHelper {

  public:
    HitHelper();
    ~HitHelper();

    // Get the track index of a track object
    size_t GetTrackIndex(const recob::Track &track, const std::vector<art::Ptr<recob::Track>> &tracklist);

    // Get the vector of art::Ptr to hit for the given track
    const std::vector<art::Ptr<recob::Hit>> GetArtPtrToHitVect(const art::FindManyP<recob::Hit> &fmht,
                                                               const size_t &trackIndex);

    // Get XYZ position for that hit
    TVector3 GetHitXYZ(const art::Ptr<recob::Hit> &hitp,
                       art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                       const std::vector<art::Ptr<recob::Track>> &tracklist,
                       const size_t &trackIndex);

    // Check if a hit has high electron contribution at a certain distance from the end point
    bool IsHitMichelLike(const art::Ptr<recob::Hit> &hitp,
                         const TVector3 &recoEndPoint,
                         art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                         const std::vector<art::Ptr<recob::Track>> &tracklist,
                         const size_t &trackIndex);


    // Set the parameters from the FHICL file
    void reconfigure(fhicl::ParameterSet const &p);

  private:
    // Fhicl file parameters for michel hits search
    double _electronEnergyFractionToCallMichelHits;
    double _maxDistanceToCallMichelHits;

    // Declare handle for particle inventory service
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Declare handle for backtracker
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // Utils
    protoana::ProtoDUNEPFParticleUtils   pfpUtil;

    // Product labels from FHICL file.
    std::string fTrackerTag, fPFParticleTag;


  };
}

#endif
