/***
  Class containing useful functions for hit on a plane.

*/
#ifndef HIT_PLANE_ALG_CXX
#define HIT_PLANE_ALG_CXX

#include "HitPlaneAlg.h"

namespace stoppingcosmicmuonselection {

  HitPlaneAlg::HitPlaneAlg(const artPtrHitVec &trackHits,
                           const size_t &start_index,
                           const size_t &planeNumber) :
                           _trackHits(trackHits),
                           _start_index(start_index),
                           _planeNumber(planeNumber) {
    _hitsOnPlane = hitHelper.GetHitsOnAPlane(_planeNumber,_trackHits);
    std::cout << "HitPlaneAlg.cxx: " << std::endl;
    std::cout << "\tSize of hits on plane before ordering: " << _hitsOnPlane.size() << std::endl;
    OrderHitVec();
    std::cout << "\tOrdered HitList Size: " << _hitsOnPlane.size()
              << std::endl << "\tWireID list size: " << _effectiveWireID.size() << std::endl;
    if (_effectiveWireID.size() != _hitsOnPlane.size())
      throw cet::exception("HitPlaneAlg.cxx") << "Hit vector and wire ID vector have different size.";
    HitSmoother();
  }

  HitPlaneAlg::~HitPlaneAlg() {

  }

  // Order hits based on their 2D (wire-time) position.
  void HitPlaneAlg::OrderHitVec() {
    std::cout << "\tOrdering hit vector..." << std::endl;
    artPtrHitVec newVector;
    newVector.reserve(_hitsOnPlane.size());
    newVector.push_back(_hitsOnPlane.at(_start_index));
    const auto &starthit = _hitsOnPlane.at(_start_index);
    _effectiveWireID.push_back(geoHelper.GetWireNumb(starthit));
    _hitsOnPlane.erase(_hitsOnPlane.begin() + _start_index);

    //double maxAllowedDistance = 50;
    //double slope_threshold = 2;

    double min_dist = DBL_MAX;
    int min_index = -1;
    std::vector<double> distances;

    while (_hitsOnPlane.size() != 0) {

      min_dist = DBL_MAX;
      min_index = -1;

      for (size_t i = 0; i < _hitsOnPlane.size(); i++) {
        // For previous hit.
        double hitPeakTime = newVector.back()->PeakTime();
        size_t wireNumb1 = geoHelper.GetWireNumb(newVector.back());
        TVector3 pt1(hitPeakTime,wireNumb1,0);
        // For current hit.
        double hitPeakTime2 = _hitsOnPlane.at(i)->PeakTime();
        size_t wireNumb2 = geoHelper.GetWireNumb(_hitsOnPlane.at(i));
        TVector3 pt2(hitPeakTime2,wireNumb2,0);
        double dist = (pt1-pt2).Mag();
        if (dist < min_dist) {
          min_index = i;
          min_dist = dist;
        }
      }

      auto const &hit = _hitsOnPlane.at(min_index);

      if (distances.size() > 3 && min_dist > (10*mean(distances))) {
        newVector.push_back(hit);
        _hitPeakTime.push_back(hit->PeakTime());
        _effectiveWireID.push_back(geoHelper.GetWireNumb(hit));
      }
      
      // std::cout << "min dist: " << min_dist << std::endl;
      // if (min_dist < maxAllowedDistance) {
      //   std::cout << "\t\tOk, adding hit. Dist < max dist." << std::endl;
      //   newVector.push_back(hit);
      //   _hitPeakTime.push_back(hit->PeakTime());
      //   _effectiveWireID.push_back(geoHelper.GetWireNumb(hit));
      // }
      // else if ((geoHelper.GetWireNumb(hit)) == _effectiveWireID.back()
      //           && min_dist < 50) {
      //   std::cout << "\t\tOk, adding hit." << std::endl;
      //   newVector.push_back(hit);
      //   _hitPeakTime.push_back(hit->PeakTime());
      //   _effectiveWireID.push_back(geoHelper.GetWireNumb(hit));
      // }
      // else if (newVector.size() > 5) {
      //   std::cout << "\t\tThe hit is too far away." << std::endl;
      //   // Calculate previous slope.
      //   auto iter = newVector.end();
      //   auto hit_2 = *(--iter);
      //   auto hit_1 = *(iter-5);
      //   double previous_slope = (hit_2->PeakTime()-hit_1->PeakTime()) / (geoHelper.GetWireNumb(hit_2)-geoHelper.GetWireNumb(hit_1));
      //   std::cout << "\t\tPrevious slope: " << previous_slope << std::endl;
      //   // Calculate next slope.
      //   double new_slope = ((hit->PeakTime()-hit_2->PeakTime()) / (geoHelper.GetWireNumb(hit)-geoHelper.GetWireNumb(hit_2)));
      //   std::cout << "\t\tCurrent slope: " << new_slope << std::endl;
      //   // Check the next hit will be in a consecutive wire
      //   bool progressive_order = false;
      //   if (geoHelper.GetWireNumb(hit_1) < geoHelper.GetWireNumb(hit_2)) {
      //     if (geoHelper.GetWireNumb(hit) > geoHelper.GetWireNumb(hit_2)) {
      //       progressive_order = true;
      //     }
      //   }
      //   if (geoHelper.GetWireNumb(hit_2) < geoHelper.GetWireNumb(hit_1)) {
      //     if (geoHelper.GetWireNumb(hit) < geoHelper.GetWireNumb(hit_2)) {
      //       progressive_order = true;
      //     }
      //   }
      //   // If the two slopes are close, then there is
      //   // probably a dead region between the point.
      //   // If so, increase the min distance by half a meter
      //   // and add the hit.
      //   if (TMath::Abs(new_slope - previous_slope) < slope_threshold &&
      //       min_dist < maxAllowedDistance + 50 &&
      //       progressive_order) {
      //     std::cout << "\t\tOk, adding hit." << std::endl;
      //     newVector.push_back(hit);
      //     _hitPeakTime.push_back(hit->PeakTime());
      //     _effectiveWireID.push_back(geoHelper.GetWireNumb(hit));
      //   }
      // } // newVector.size() > 5

      _hitsOnPlane.erase(_hitsOnPlane.begin() + min_index);
    }
    _areHitOrdered = true;
    std::swap(_hitsOnPlane, newVector);
    std::cout << "\tHit vector ordered." << std::endl;
    return;
  }

  // Smooth hits.
  void HitPlaneAlg::HitSmoother() {
    if (!_areHitOrdered)
      OrderHitVec();
    std::cout << "\tSmoothing hits..." << std::endl;
    artPtrHitVec newVector;
    std::vector<double> newVector_wire;
    std::vector<double> meanVec;
    std::vector<double> wireVec;
    if (_effectiveWireID.size()<=2) return;
    for (const auto &bunch : get_neighbors(_effectiveWireID, 2)) {
      for (const auto &wire : bunch) {
        wireVec.push_back(wire);
      }
      meanVec.push_back(mean(wireVec));
      wireVec.clear();
    }
    newVector.push_back(_hitsOnPlane.at(0));
    newVector.push_back(_hitsOnPlane.at(1));
    newVector_wire.push_back(_effectiveWireID.at(0));
    newVector_wire.push_back(_effectiveWireID.at(1));
    std::cout << "\t\tGot wire mean..." << std::endl;
    std::cout << "\t\tWireid size: " << _effectiveWireID.size() << std::endl;
    std::cout << "\t\thitsOnPlane size: " << _hitsOnPlane.size() << std::endl;
    std::cout << "\t\tmeanVec size: " << meanVec.size() << std::endl;

    for (size_t i = 2; i < meanVec.size()-1; i++) {
      if (std::abs(meanVec.at(i-1) - meanVec.at(i)) < 1        &&
          _effectiveWireID.at(i-1) !=  _effectiveWireID.at(i)  &&
          std::abs(meanVec.at(i)   - meanVec.at(i+1)) < 1      &&
          _effectiveWireID.at(i) !=  _effectiveWireID.at(i+1) ) {
        //std::cout << "\t\tIn the if" << std::endl;
        if (_hitsOnPlane.at(i)->Integral() > _hitsOnPlane.at(i+1)->Integral()) {
          newVector.push_back(_hitsOnPlane.at(i));
          newVector_wire.push_back(_effectiveWireID.at(i));
        }
        else {
          newVector.push_back(_hitsOnPlane.at(i+1));
          newVector_wire.push_back(_effectiveWireID.at(i+1));
        }

        i++;
      }
      else {
        newVector.push_back(_hitsOnPlane.at(i));
        newVector_wire.push_back(_effectiveWireID.at(i));
      }
    }
    std::cout << "\tSmoothing ok." << std::endl;
    newVector.push_back(_hitsOnPlane.at(_hitsOnPlane.size()-1));
    newVector_wire.push_back(_effectiveWireID.at(_hitsOnPlane.size()-1));
    std::swap(newVector, _hitsOnPlane);
    std::swap(newVector_wire, _effectiveWireID);
  }

  // Get the ordered hit vector.
  const artPtrHitVec HitPlaneAlg::GetOrderedHitVec() {
    if (!_areHitOrdered)
      OrderHitVec();
    return _hitsOnPlane;
  }

  // Get the ordered wire number.
  const std::vector<double> HitPlaneAlg::GetOrderedWireNumb() {
    if (!_areHitOrdered)
      OrderHitVec();
    return _effectiveWireID;
  }

  // Work out the vector of ordered hit charge.
  const std::vector<double> HitPlaneAlg::GetOrderedQ() {
    if (!_areHitOrdered)
      OrderHitVec();
    std::vector<double> Qs;
    //std::cout << "Calculating hit Qs..." << std::endl;
    for (size_t i = 0; i < _hitsOnPlane.size()-1; i++)
      Qs.push_back(_hitsOnPlane[i]->Integral());
    return Qs;
  }

  // Work out the vector of ordered dQds.
  const std::vector<double> HitPlaneAlg::GetOrderedDqds() {
    if (!_areHitOrdered)
      OrderHitVec();
    std::vector<double> dQds;
    double ds = 1.;
    //std::cout << "Calculating dQds..." << std::endl;
    for (size_t i = 0; i < _hitsOnPlane.size()-1; i++) {
      auto const &hit = _hitsOnPlane[i];
      auto const &nextHit = _hitsOnPlane[i+1];
      double XThisPoint = detprop->ConvertTicksToX(hit->PeakTime(),hit->WireID().Plane,hit->WireID().TPC,hit->WireID().Cryostat);
      double XNextPoint = detprop->ConvertTicksToX(nextHit->PeakTime(),nextHit->WireID().Plane,nextHit->WireID().TPC,nextHit->WireID().Cryostat);

      TVector3 thisPoint(_effectiveWireID[i]*geoHelper.GetWirePitch(_planeNumber), XThisPoint, 0);
      TVector3 nextPoint(_effectiveWireID[i+1]*geoHelper.GetWirePitch(_planeNumber), XNextPoint, 0);

      ds = (thisPoint - nextPoint).Mag();
      dQds.push_back(hit->Integral() / ds);
    }
    // Need final point.
    dQds.push_back(_hitsOnPlane.at(_hitsOnPlane.size()-1)->Integral() / ds);
    //std::cout << "Dqds vector for this track calculated." << std::endl;
    return dQds;
  }

  // Define smoother.
  const std::vector<double> HitPlaneAlg::Smoother(const std::vector<double> &object, const size_t &Nneighbors) {
    std::vector<double> result;
    result.reserve(object.size());
    for (const auto &neighbors : get_neighbors(object,Nneighbors)) {
      double median = get_smooth_trunc_median(neighbors);
      result.push_back(median);
    }
    return result;
  }

  // Calculate local linearity.
  const std::vector<double> HitPlaneAlg::CalculateLocalLinearity(const size_t &Nneighbors) {
    std::vector<double> linearity;
    std::vector<double> time, wire;
    for (const auto &hits : get_neighbors(_hitsOnPlane,Nneighbors)) {
      for (const auto &hit : hits) {
        time.push_back(hit->PeakTime());
        wire.push_back(geoHelper.GetWireNumb(hit));
      }
      double covariance = cov(time,wire);
      double stdevTime = stdev(time);
      double stdevWire = stdev(wire);
      double N = time.size();
      double lin = TMath::Abs(covariance) / N / (stdevTime*stdevWire);
      if (std::isnan(lin)) lin = 0.0;
      linearity.push_back(lin);
      time.clear();
      wire.clear();
    }
    _isLinearityCalculated = true;
    return linearity;
  }

} // end of namespace stoppingcosmicmuonselection

#endif
