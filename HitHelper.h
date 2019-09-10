/***
  Class containing useful functions for hit.

*/
#ifndef HIT_HELPER_H
#define HIT_HELPER_H

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

    // Get track index in track vector for that event.
    size_t GetTrackIndex(const recob::Track &track, const std::vector<art::Ptr<recob::Track>> &tracklist);

    // Get XYZ position for that hit
    TVector3 GetHitXYZ(const art::Ptr<recob::Hit> &hitp,
                       art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                       std::vector<art::Ptr<recob::Track>> &tracklist,
                       size_t &trackIndex);

    //
    bool IsTrackWithMichelHits();

    // Set the parameters from the FHICL file
    void reconfigure(fhicl::ParameterSet const &p);

  private:

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
