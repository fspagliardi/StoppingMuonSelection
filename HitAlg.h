/***
  Class containing useful algorithms for hit.

*/
#ifndef HIT_ALG_H
#define HIT_ALG_H

#include "TVector3.h"
#include "TMath.h"

#include "DataTypes.h"
#include "HitHelper.h"

namespace stoppingcosmicmuonselection {

  class HitAlg {

  public:
    HitAlg(const std::vector<art::Ptr<recob::Hit>> &trackHits,
           art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
           const std::vector<art::Ptr<recob::Track>> &tracklist,
           const size_t &trackIndex);
    ~HitAlg();

    // Get hit list on a given plane.
    std::vector<art::Ptr<recob::Hit>> GetHitsOnAPlane(const size_t &planeNumb);

    // Get the closest hit to a given point on the given track.
    art::Ptr<recob::Hit> GetClosestHitToPoint(const TVector3 &point, const size_t &planeNumb);

    // Get index of the closest hit to a given point on the given track.
    const size_t GetIndexClosestHitToPoint(const TVector3 &point, const size_t &planeNumb);

    // Order hits based on their 2D (wire-time) position.
    const std::vector<art::Ptr<recob::Hit>> Get2DOrderedHitArray(const TVector3 &startPoint,
                                                                         const size_t &planeNumb);

  private:

    std::vector<art::Ptr<recob::Hit>> _hitsOnPlane;

    // To be initialised.
    const std::vector<art::Ptr<recob::Hit>> &_trackHits;
    art::FindManyP<recob::Hit,recob::TrackHitMeta> &_fmthm;
    const std::vector<art::Ptr<recob::Track>> &_tracklist;
    const size_t &_trackIndex;

    // Declare the hit helper.
    HitHelper hitHelper;

  };
}

#endif
