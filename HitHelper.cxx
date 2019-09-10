/***
  Class containing useful functions for hit.

*/
#ifndef HIT_HELPER_CXX
#define HIT_HELPER_CXX

#include "HitHelper.h"

namespace stoppingcosmicmuonselection {

  HitHelper::HitHelper() {

  }

  HitHelper::~HitHelper() {

  }

  // Get the track index of a track object
  size_t HitHelper::GetTrackIndex(const recob::Track &track,
                                  const std::vector<art::Ptr<recob::Track>> &tracklist) {
    if (tracklist.size() == 0) return INV_INT;
    for (size_t trackIter = 0; trackIter < tracklist.size(); trackIter++) {
      if (tracklist[trackIter]->ID() == track.ID())
        return trackIter;
    }
    return INV_INT;
  }

  // Get the vector of art::Ptr to hit for the given track
  const std::vector<art::Ptr<recob::Hit>> HitHelper::GetArtPtrToHitVect(const art::FindManyP<recob::Hit> &fmht,
                                                                        const size_t &trackIndex) {
    return fmht.at(trackIndex);
  }

  // Get XYZ position for that hit
  TVector3 HitHelper::GetHitXYZ(const art::Ptr<recob::Hit> &hitp,
                                art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                const std::vector<art::Ptr<recob::Track>> &tracklist,
                                const size_t &trackIndex) {

    TVector3 hitLoc(INV_DBL,INV_DBL,INV_DBL);
    if (!fmthm.isValid()) return hitLoc;

    const std::vector<art::Ptr<recob::Hit>> vhit = fmthm.at(trackIndex);
    const std::vector<const recob::TrackHitMeta*> vmeta = fmthm.data(trackIndex);
    // iterate on meta data
    for (size_t ii=0;ii<vhit.size();++ii) {
      if (vhit[ii].key() != hitp.key())
        continue;
      if (vmeta[ii]->Index() == std::numeric_limits<int>::max()) {
        continue;
      }
      if (vmeta[ii]->Index()>=tracklist[trackIndex]->NumberTrajectoryPoints()){
        throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[trackIndex]->NumberTrajectoryPoints()<<" for track index "<<trackIndex<<". Something is wrong";
      }
      if (!tracklist[trackIndex]->HasValidPoint(vmeta[ii]->Index())){
        continue;
      }
      hitLoc = tracklist[trackIndex]->LocationAtPoint<TVector3>(vmeta[ii]->Index());
    } // iteration on metadata

    return hitLoc;
  }

  // Check if a hit has high electron contribution at a certain distance from the end point
  bool HitHelper::IsHitMichelLike(const art::Ptr<recob::Hit> &hitp,
                                  const TVector3 &recoEndPoint,
                                  art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                  const std::vector<art::Ptr<recob::Track>> &tracklist,
                                  const size_t &trackIndex) {
    for (const sim::TrackIDE &tIDE : bt_serv->HitToTrackIDEs(*hitp)) {
      // check if there is a contribution from an electron.
      if (TMath::Abs(pi_serv->TrackIdToParticle_P(tIDE.trackID)->PdgCode())==11) {
        TVector3 hitLoc = GetHitXYZ(hitp,fmthm,tracklist,trackIndex);
        if (hitLoc == TVector3(INV_DBL,INV_DBL,INV_DBL)) continue;
        if (tIDE.energyFrac>_electronEnergyFractionToCallMichelHits
            && (hitLoc-recoEndPoint).Mag()<_maxDistanceToCallMichelHits)
          return true;
      }
    }
    return false;
  }

  // Set the parameters from the FHICL file
  void HitHelper::reconfigure(fhicl::ParameterSet const &p) {
    fTrackerTag = p.get<std::string>("TrackerTag");
    fPFParticleTag = p.get<std::string>("PFParticleTag");
    _electronEnergyFractionToCallMichelHits = p.get<double>("electronEnergyFractionToCallMichelHits", 0.7);
    _maxDistanceToCallMichelHits = p.get<double>("maxDistanceToCallMichelHits", 15);
  }

} // end of namespace stoppingcosmicmuonselection

#endif
