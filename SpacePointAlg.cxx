/***
  Class containing useful functions for not fitted space points.

*/
#ifndef SPACEPOINT_ALG_CXX
#define SPACEPOINT_ALG_CXX

#include "SpacePointAlg.h"

namespace stoppingcosmicmuonselection {

SpacePointAlg::SpacePointAlg() {
}

SpacePointAlg::~SpacePointAlg() {

}

// Is valid for this track
bool SpacePointAlg::IsValid() {
  return _isValid;
}

// Read parameters from FHICL file
void SpacePointAlg::reconfigure(fhicl::ParameterSet const &p) {
  _cilinderAxis = p.get<double>("cilinderAxis", 50.);
  _cilinderRadius = p.get<double>("cilinderRadius", 5.);
  _minNumberSpacePoints = p.get<int>("minNumberSpacePoints", 10);
}

// Check if the track correctly fit the space points around the end points
bool SpacePointAlg::IsGoodTrack(const recob::Track &track,
                                const std::vector<recob::SpacePoint> &spacePoints,
                                const trackProperties &trackProp) {
  tP = trackProp;
  // Check if the track is valid
  bool isTrackValid = false;
  size_t nValidPointsCounter(0);
  for (size_t j = 0; j < track.NumberTrajectoryPoints(); j++)   {
    if (track.HasValidPoint(j))
      nValidPointsCounter++;
  }
  if (nValidPointsCounter > 2)
    isTrackValid = true;

  if (isTrackValid && track.Length()>100.)
    _isValid = true;

  TVector3 posLastValidPoint, pos20cmLastValidPoint, posFirstValidPoint, pos20cmFirstValidPoint;
  if (_isValid) {
    // Initialise points
    posLastValidPoint = track.LocationAtPoint<TVector3>(track.LastValidPoint());
    posFirstValidPoint = track.LocationAtPoint<TVector3>(track.FirstValidPoint());
    bool found1 = false, found2 = false;

    for (size_t ii = 0; ii <= track.NumberTrajectoryPoints(); ii++)   {
      if (!track.HasValidPoint(ii)) continue;
      pos20cmLastValidPoint = track.LocationAtPoint<TVector3>(ii);
      if (TMath::Abs((pos20cmLastValidPoint-posLastValidPoint).Mag()-20.) < 5) {
        found1 = true;
        break;
      }
    }
    for (size_t ii = 0; ii <= track.NumberTrajectoryPoints(); ii++)   {
      if (!track.HasValidPoint(ii)) continue;
      pos20cmFirstValidPoint = track.LocationAtPoint<TVector3>(ii);
      if (TMath::Abs((pos20cmFirstValidPoint-posFirstValidPoint).Mag()-20.) < 5) {
        found2 = true;
        break;
      }
    }
    //std::cout << found1 << " " << found2 << std::endl;
    if (!found1 || !found2) {
      pos20cmLastValidPoint = posFirstValidPoint;
      pos20cmFirstValidPoint = posLastValidPoint;
    }

    bool badTrack = true;
    if (tP.isAnodeCrosserMine)
      badTrack = (IsTrackNotFittingSpacePoints(posLastValidPoint,pos20cmLastValidPoint,spacePoints,_bottom)
           || IsTrackNotFittingSpacePoints(posFirstValidPoint,pos20cmFirstValidPoint,spacePoints,_top));
    else if (tP.isAnodeCrosserPandora || tP.isCathodeCrosser)
      badTrack = IsTrackNotFittingSpacePoints(posLastValidPoint,pos20cmLastValidPoint,spacePoints,_bottom);

    return !badTrack;
  } // if track is valid
  else
    return false;
}

// Given a point and a line find the projection of that point on the line in 2D
TVector3 SpacePointAlg::FindFoot(double *coeffLine, const double &sp_Y, const double &sp_Z)  {
  TVector3 foot;
  double zz = (sp_Y + (sp_Z / coeffLine[1]) - coeffLine[0]) / ( ((coeffLine[1]*coeffLine[1])+1) / coeffLine[1] );
  double yy = (coeffLine[1] * zz) + coeffLine[0];
  foot.SetXYZ(INV_DBL,yy,zz);
  return foot;
}

// Check if the track is missing some space points for one end
void SpacePointAlg::FillLineCoeff(TVector3 &posLastValidPoint,
                                  TVector3 &pos20cmLastValidPoint,
                                  double *coeffLineYZ,
                                  double *coeffLineXZ) {
  coeffLineYZ[1] = (posLastValidPoint.Y()-pos20cmLastValidPoint.Y()) / (posLastValidPoint.Z()-pos20cmLastValidPoint.Z());
  coeffLineYZ[0] = pos20cmLastValidPoint.Y() - (pos20cmLastValidPoint.Z()*(coeffLineYZ[1]));
  coeffLineXZ[1] = (posLastValidPoint.X()-pos20cmLastValidPoint.X()) / (posLastValidPoint.Z()-pos20cmLastValidPoint.Z());
  coeffLineXZ[0] = pos20cmLastValidPoint.X() - (pos20cmLastValidPoint.Z()*(coeffLineXZ[1]));
}

// Check if the track is missing some space points for one end
bool SpacePointAlg::IsTrackNotFittingSpacePoints(TVector3 &posExtremeValidPoint,
                                                    TVector3 &pos20cmValidPoint,
                                                    const std::vector<recob::SpacePoint> &spacePoints,
                                                    const std::string &whichEnd) {
  //std::cout << "Working with option ---> " << whichEnd << std::endl;
  double coeffLineYZ[2] = {INV_DBL,INV_DBL}, coeffLineXZ[2] = {INV_DBL,INV_DBL};
  FillLineCoeff(posExtremeValidPoint,pos20cmValidPoint,coeffLineYZ,coeffLineXZ);
  //std::cout << coeffLineYZ[0] << " " << coeffLineYZ[1] << " " << coeffLineXZ[0] << " " << coeffLineXZ[1] << std::endl;
  //std::cout << posExtremeValidPoint.Y() << " " << posExtremeValidPoint.Z() << std::endl;
  //std::cout << pos20cmValidPoint.Y() << " " << pos20cmValidPoint.Z() << std::endl;
  size_t spCounter = 0; // count SP

  // Iterates on space points
  for (auto const & sp : spacePoints) {

    // Now look at stuff in YZ plane
    TVector3 footYZ = FindFoot(coeffLineYZ,sp.XYZ()[1],sp.XYZ()[2]);
    // Now look at stuff in XZ plane
    TVector3 footXZ = FindFoot(coeffLineXZ,sp.XYZ()[0],sp.XYZ()[2]);

    // Check if the space point is within the designed geometry. NB: The coordinated for the vector foot**
    // are always such that Y(Z) for every plane
    double distanceFootSpYZ = TMath::Sqrt(TMath::Power(footYZ.Y()-sp.XYZ()[1],2) + TMath::Power(footYZ.Z()-sp.XYZ()[2],2));
    double distanceFootEndYZ = TMath::Sqrt(TMath::Power(footYZ.Y()-posExtremeValidPoint.Y(),2) + TMath::Power(footYZ.Z()-posExtremeValidPoint.Z(),2));
    double distanceFootSpXZ = TMath::Sqrt(TMath::Power(footXZ.Y()-sp.XYZ()[0],2) + TMath::Power(footXZ.Z()-sp.XYZ()[2],2));
    double distanceFootEndXZ = TMath::Sqrt(TMath::Power(footXZ.Y()-posExtremeValidPoint.X(),2) + TMath::Power(footXZ.Z()-posExtremeValidPoint.Z(),2));

    if (whichEnd == "bottom") { // option bottom includes cathode-crossing tracks as well
      // If track is T0 tagged the space points not fitted are not aligned in the XY place because they have
      // the wrong X while the T0-tagged track has been shifted.
      if (tP.isAnodeCrosserPandora || tP.isCathodeCrosser) {
        if (distanceFootSpYZ<_cilinderRadius && distanceFootEndYZ<_cilinderAxis && posExtremeValidPoint.Y()>sp.XYZ()[1]) {
          spCounter++;
        }
      }
      else if (tP.isAnodeCrosserMine) {
        if (distanceFootSpYZ<_cilinderRadius && distanceFootEndYZ<_cilinderAxis && distanceFootSpXZ<_cilinderRadius && distanceFootEndXZ<_cilinderAxis && posExtremeValidPoint.Y()>sp.XYZ()[1]) {
          spCounter++;
          //std::cout << "X: " << sp.XYZ()[0] << " Y: " << sp.XYZ()[1] << " Z: " << sp.XYZ()[2] << std::endl;
        }
      }

    }
    else if (whichEnd == "top") {
      if (tP.isAnodeCrosserMine) {
        if (distanceFootSpYZ<_cilinderRadius && distanceFootEndYZ<_cilinderAxis && distanceFootSpXZ<_cilinderRadius && distanceFootEndXZ<_cilinderAxis && posExtremeValidPoint.Y()<sp.XYZ()[1]) {
          spCounter++;
          //std::cout << "X: " << sp.XYZ()[0] << " Y: " << sp.XYZ()[1] << " Z: " << sp.XYZ()[2] << std::endl;
        }
      }
    }
    else {
      //std::cout << "ERROR: wrong option for function IsTrackNotFittingSpacePoints" << std::endl;
      break;
    }

  } // end loop over space points

  //std::cout << "Number of SP in geometry: " << spCounter << std::endl;
  if (spCounter > _minNumberSpacePoints)
    return true;
  else
    return false;
}

} // end of namespace stoppingcosmicmuonselection

#endif
