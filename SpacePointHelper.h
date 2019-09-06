/***
  Class containing useful functions for not fitted space points.

*/
#ifndef SPACEPOINT_HELPER_H
#define SPACEPOINT_HELPER_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "TVector3.h"
#include "TMath.h"

namespace stoppingcosmicmuonselection {

  class SpacePointHelper {

  public:
    SpacePointHelper(const double &axis, const double &radius, const double &minNumber, const double &T0);
    ~SpacePointHelper();

    // Is valid for this track
    bool IsValid();

    // Check if the track correctly fit the space points around the end points
    bool IsGoodTrack(const recob::Track &track, std::vector<recob::SpacePoint> &spacePoints);

  private:
    // Strings for track orientation
    const std::string _bottom = "bottom";
    const std::string _top = "top";
    double _T0;
    bool _isValid = false;
    // Define parameters for cilinder space points
    double _cilinderAxis; // cm
    double _cilinderRadius; // cm
    size_t _minNumberSpacePoints;
    const double INV_DBL = -9999999;
    // Given a point and a line find the projection of that point on the line in 2D
    TVector3 FindFoot(double *coeffLine, const double &sp_Y, const double &sp_Z);
    // Fill the coefficients for the line interpolating the start and end of the track
    void FillLineCoeff(TVector3 &posLastValidPoint,
                       TVector3 &pos20cmLastValidPoint,
                       double *coeffLineYZ,
                       double *coeffLineXZ);
    // Check if the track is missing some space points for one end
    bool IsTrackNotFittingSpacePoints(TVector3 &posExtremeValidPoint,
                                      TVector3 &pos20cmValidPoint,
                                      const std::vector<recob::SpacePoint> &spacePoints,
                                      const std::string &whichEnd);

  };
} // namespace

#endif
