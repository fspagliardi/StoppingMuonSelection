/***
  Class containing useful functions for hit.

*/
#ifndef HIT_ALG_CXX
#define HIT_ALG_CXX

#include "HitAlg.h"

namespace stoppingcosmicmuonselection {

  HitAlg::HitAlg(const std::vector<art::Ptr<recob::Hit>> &trackHits,
                 art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                const std::vector<art::Ptr<recob::Track>> &tracklist,
                const size_t &trackIndex) {

  }

  HitAlg::~HitAlg() {

  }

  // Get hit list on a given plane.
  std::vector<art::Ptr<recob::Hit>> HitAlg::GetHitsOnAPlane(const size_t &planeNumb) {
    std::vector<art::Ptr<recob::Hit>> hitsOnPlane;
    for (auto const &hitp : _trackHits) {
      if (!hitp->WireID().isValid) continue;
      if (hitp->WireID().Plane != planeNumb) continue;
      hitsOnPlane.push_back(hitp);
    }
    return hitsOnPlane;
  }

  // Get the closest hit to a given point on the given track.
  art::Ptr<recob::Hit> HitAlg::GetClosestHitToPoint(const TVector3 &point, const size_t &planeNumb) {
    size_t hitIndex = GetIndexClosestHitToPoint(point, planeNumb);
    return _trackHits.at(hitIndex);
  }

  // Get index of the closest hit to a given point on the given track.
  const size_t HitAlg::GetIndexClosestHitToPoint(const TVector3 &point, const size_t &planeNumb) {
    std::vector<double> distanceVec;
    for (auto const &hitp : GetHitsOnAPlane(planeNumb)) {
      const TVector3 &hitLoc = hitHelper.GetHitXYZ(hitp,_fmthm,_tracklist,_trackIndex);
      distanceVec.push_back((hitLoc - point).Mag());
    }
    auto it_min_element = std::min_element(distanceVec.begin(),distanceVec.end());
    size_t hitIndex = it_min_element - distanceVec.begin();
    return hitIndex;
  }

  // Order hits based on their 2D (wire-time) position.
  const std::vector<art::Ptr<recob::Hit>> HitAlg::Get2DOrderedHitArray(const TVector3 &startPoint,
                                                                       const size_t &planeNumb) {
    std::vector<art::Ptr<recob::Hit>> hitVectorOrdered;
    std::vector<art::Ptr<recob::Hit>> hitsOnPlane = GetHitsOnAPlane(planeNumb);
    hitVectorOrdered.reserve(hitsOnPlane.size());
    // Push back starting hit.
    art::Ptr<recob::Hit> startHit = GetClosestHitToPoint(startPoint,planeNumb);
    hitVectorOrdered.push_back(startHit);
    //hitsOnPlane.erase(hitsOnPlane.at(GetIndexClosestHitToPoint(startPoint,planeNumb)));

    double min_dist = DBL_MAX;
    int min_index = -1;

    while (hitsOnPlane.size() != 0) {

      min_dist = DBL_MAX;
      min_index = -1;

      for (size_t i = 0; i < hitsOnPlane.size(); i++) {

      }
    }
    return hitVectorOrdered;
  }

} // end of namespace stoppingcosmicmuonselection

#endif
