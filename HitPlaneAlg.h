/***
  Class containing useful algorithms for hit on a plane.

*/
#ifndef HIT_PLANE_ALG_H
#define HIT_PLANE_ALG_H

#include "TVector3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "DataTypes.h"
#include "HitHelper.h"
#include "GeometryHelper.h"
#include "Tools.h"

namespace stoppingcosmicmuonselection {

  class HitPlaneAlg {

  public:
    HitPlaneAlg(const artPtrHitVec &trackHits, const size_t &start_index, const size_t &planeNumber);
    ~HitPlaneAlg();

    // Order hits based on their 2D (wire-time) position.
    void OrderHitVec();

    // Smooth hits.
    void HitSmoother();

    // Get the ordered hit vector.
    const artPtrHitVec GetOrderedHitVec();

    // Get the ordered wire number.
    const std::vector<double> GetOrderedWireNumb();

    // Work out the vector of ordered hit charge.
    const std::vector<double> GetOrderedQ();

    // Work out the vector of ordered dQds.
    const std::vector<double> GetOrderedDqds();

    // Define smoother.
    const std::vector<double> Smoother(const std::vector<double> &object, const size_t &Nneighbors);

    // Calculate local linearity.
    const std::vector<double> CalculateLocalLinearity(const size_t &Nneighbors);

  private:
    const artPtrHitVec &_trackHits;
    const size_t &_start_index;
    const size_t &_planeNumber;

    artPtrHitVec _hitsOnPlane;
    std::vector<double> _hitPeakTime;
    std::vector<double> _effectiveWireID;

    bool _areHitOrdered = false;
    bool _isLinearityCalculated = false;

    // Helpers.
    GeometryHelper geoHelper;
    HitHelper      hitHelper;

    // Declare handle for detector properties
    const detinfo::DetectorProperties *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  };
}

#endif
