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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "TVector3.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TGraph2D.h"

#include "DataTypes.h"
#include "GeometryHelper.h"

namespace stoppingcosmicmuonselection {

  class HitHelper {

  public:
    HitHelper();
    ~HitHelper();

    // Get the track index of a track object
    size_t GetTrackIndex(const recob::Track &track, const std::vector<art::Ptr<recob::Track>> &tracklist);

    // Get the vector of art::Ptr to hit for the given track
    const artPtrHitVec GetArtPtrToHitVect(const art::FindManyP<recob::Hit> &fmht,
                                          const size_t &trackIndex);

    // Get hit list on a given plane.
    artPtrHitVec GetHitsOnAPlane(const size_t &planeNumb,
                                 const artPtrHitVec &allHits);

    // Get XYZ position for that hit
    TVector3 GetHitXYZ(const art::Ptr<recob::Hit> &hitp,
                       art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                       const std::vector<art::Ptr<recob::Track>> &tracklist,
                       const size_t &trackIndex);

    // Get index of the closest hit to a given point on the given track.
    const size_t GetIndexClosestHitToPoint(const TVector3 &point,
                                           const artPtrHitVec &hits,
                                           art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                           const std::vector<art::Ptr<recob::Track>> &tracklist,
                                           const size_t &trackIndex);

    // Get the closest hit to a given point on the given track.
    art::Ptr<recob::Hit> GetClosestHitToPoint(const TVector3 &point,
                                              const artPtrHitVec &hits,
                                              art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                              const std::vector<art::Ptr<recob::Track>> &tracklist,
                                              const size_t &trackIndex);

    // Check if a hit has high electron contribution at a certain distance from the end point
    bool IsHitMichelLike(const art::Ptr<recob::Hit> &hitp,
                         const TVector3 &recoEndPoint,
                         art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                         const std::vector<art::Ptr<recob::Track>> &tracklist,
                         const size_t &trackIndex);

    // Get subvector of michel-like hits.
    artPtrHitVec GetMichelLikeHits(const artPtrHitVec &hits,
                                   const TVector3 &recoEndPoint,
                                   art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                   const std::vector<art::Ptr<recob::Track>> &tracklist,
                                   const size_t &trackIndex);

    // Get subvector of michel-like hits.
    artPtrHitVec GetMuonLikeHits(const artPtrHitVec &hits,
                                 const TVector3 &recoEndPoint,
                                 art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                 const std::vector<art::Ptr<recob::Track>> &tracklist,
                                 const size_t &trackIndex);

    // Fill the TGraph2D for the images.
    void FillTrackGraph2D(TGraph2D *graph,
                          const artPtrHitVec &trackHits,
                          const TVector3 &recoEndPoint,
                          const size_t &planeNumber,
                          const double &t0);

    // Get a TProfile2D filled with hit peak times and wire number
    void FillTrackHitPicture(TProfile2D* image,
                             const artPtrHitVec &trackHits,
                             const TVector3 &recoEndPoint,
                             const size_t &planeNumber);

    // Initialise the image for a series of hit for a given plane
    void InitHitImageHisto(TProfile2D *&image, const size_t &planeNumber, const std::string &name);

    // Check if vector of hits contain hits on the cryostat side.
    bool AreThereHitsOnCryoSide(const std::vector<const recob::Hit*> &hits);

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

    // Declare handle for detector properties
    const detinfo::DetectorProperties *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Geometry helper.
    GeometryHelper geoHelper;

  };
}

#endif
