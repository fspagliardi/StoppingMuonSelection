/***
  Class containing useful functions for hit on a plane.

*/
#ifndef HIT_PLANE_ALG_CXX
#define HIT_PLANE_ALG_CXX

#include "HitPlaneAlg.h"

namespace stoppingcosmicmuonselection {

  HitPlaneAlg::HitPlaneAlg(artPtrHitVec &hitsOnPlane,
                           const size_t &start_index,
                           const size_t &planeNumber) :
                           _hitsOnPlane(hitsOnPlane),
                           _start_index(start_index),
                           _planeNumber(planeNumber) {
  }

  HitPlaneAlg::~HitPlaneAlg() {

  }

  // Order hits based on their 2D (wire-time) position.
  const artPtrHitVec HitPlaneAlg::GetOrderedHitArray() {
    artPtrHitVec newVector;
    newVector.reserve(_hitsOnPlane.size());
    newVector.push_back(_hitsOnPlane.at(_start_index));
    _hitsOnPlane.erase(_hitsOnPlane.begin() + _start_index);

    double min_dist = DBL_MAX;
    int min_index = -1;

    while (_hitsOnPlane.size() != 0) {

      min_dist = DBL_MAX;
      min_index = -1;

      for (size_t i = 0; i < _hitsOnPlane.size(); i++) {
        // For previous hit.
        unsigned int hit_tpcid = newVector.back()->WireID().TPC;
        double hitPeakTime = newVector.back()->PeakTime();
        unsigned int wireID = newVector.back()->WireID().Wire;
        size_t wireOffset = geoHelper.GetWireOffset(hit_tpcid, _planeNumber);
        TVector3 pt1(hitPeakTime,wireID+wireOffset,0);
        // For current hit.
        unsigned int hit_tpcid2 = _hitsOnPlane.at(i)->WireID().TPC;
        double hitPeakTime2 = _hitsOnPlane.at(i)->PeakTime();
        unsigned int wireID2 = _hitsOnPlane.at(i)->WireID().Wire;
        size_t wireOffset2 = geoHelper.GetWireOffset(hit_tpcid2, _planeNumber);
        TVector3 pt2(hitPeakTime2,wireID2+wireOffset2,0);
        double dist = (pt1-pt2).Mag();
        if (dist < min_dist) {
          min_index = i;
          min_dist = dist;
        }
      }
      newVector.push_back(_hitsOnPlane.at(min_index));
      _hitsOnPlane.erase(_hitsOnPlane.begin() + min_index);
    }
    _isOrdered = true;
    _hitsOnPlane = newVector;
    return _hitsOnPlane;
  }

} // end of namespace stoppingcosmicmuonselection

#endif
