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

  //
  size_t HitHelper::GetTrackIndex(const recob::Track &track,
                                  const std::vector<art::Ptr<recob::Track>> &tracklist) {
    if (tracklist.size() == 0) return INV_INT;
    for (size_t trackIter = 0; trackIter < tracklist.size(); trackIter++) {
      if (tracklist[trackIter]->ID() == track.ID())
        return trackIter;
    }
    return INV_INT;
  }

  // Get XYZ position for that hit
  TVector3 GetHitXYZ(const art::Ptr<recob::Hit> &hitp,
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

  // Set the parameters from the FHICL file
  void HitHelper::reconfigure(fhicl::ParameterSet const &p) {
    fTrackerTag = p.get<std::string>("TrackerTag");
    fPFParticleTag = p.get<std::string>("PFParticleTag");
  }

} // end of namespace stoppingcosmicmuonselection

#endif
