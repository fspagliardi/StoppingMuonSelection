/***
  Class containing useful functions for hit on a plane.

*/
#ifndef HIT_PLANE_ALG_CXX
#define HIT_PLANE_ALG_CXX

#include "HitPlaneAlg.h"

namespace stoppingcosmicmuonselection {

  HitPlaneAlg::HitPlaneAlg(artPtrHitVec &hitsOnPlane,
                           const size_t start_index) :
                           _hitsOnPlane(hitsOnPlane),
                           _start_index(start_index){
  }

  HitPlaneAlg::~HitPlaneAlg() {

  }

  // Order hits based on their 2D (wire-time) position.
  // const artPtrHitVec HitPlaneAlg::GetOrderedHitArray() {
  //   artPtrHitVec hitVectorOrdered;
  //   hitVectorOrdered.reserve(_hitsOnPlane.size());
  //   //hitsOnPlane.erase(hitsOnPlane.at(GetIndexClosestHitToPoint(startPoint,planeNumb)));
  //
  //   double min_dist = DBL_MAX;
  //   int min_index = -1;
  //
  //   while (_hitsOnPlane.size() != 0) {
  //
  //     min_dist = DBL_MAX;
  //     min_index = -1;
  //
  //     for (size_t i = 0; i < _hitsOnPlane.size(); i++) {
  //
  //     }
  //   }
  //   return hitVectorOrdered;
  // }

} // end of namespace stoppingcosmicmuonselection

#endif
