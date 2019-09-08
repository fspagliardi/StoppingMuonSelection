/***
  Class containing useful functions for geometry.

*/
#ifndef GEOMETRY_HELPER_CXX
#define GEOMETRY_HELPER_CXX

#include "GeometryHelper.h"

namespace stoppingcosmicmuonselection {

  GeometryHelper::GeometryHelper() {

  }

  GeometryHelper::~GeometryHelper() {

  }

  // Initialise active volume
  void GeometryHelper::InitActiveVolumeBounds() {
    _activeBounds[0] = _activeBounds[2] = _activeBounds[4] = DBL_MAX;
    _activeBounds[1] = _activeBounds[3] = _activeBounds[5] = -DBL_MAX;
    const geo::GeometryCore *geom = lar::providerFrom<geo::Geometry>();
    double abs_X_collection = 0;
    for (geo::TPCGeo const& TPC: geom->IterateTPCs())  {
      double origin[3] = {0.};
      double center[3] = {0.};
      TPC.LocalToWorld(origin, center); // had to modify CMakeLists.txt to make this work
      double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length()};
      if( center[0] - tpcDim[0] < _activeBounds[0]) _activeBounds[0] = center[0] - tpcDim[0];
      if( center[0] - tpcDim[0] > _activeBounds[1]) _activeBounds[1] = center[0] + tpcDim[0];
      if( center[1] - tpcDim[1] < _activeBounds[2]) _activeBounds[2] = center[1] - tpcDim[1];
      if( center[1] - tpcDim[1] > _activeBounds[3]) _activeBounds[3] = center[1] + tpcDim[1];
      if( center[2] - tpcDim[2] < _activeBounds[4]) _activeBounds[4] = center[2] - tpcDim[2];
      if( center[2] - tpcDim[2] > _activeBounds[5]) _activeBounds[5] = center[2] + tpcDim[2];
      //check coordinates of collection plane
      geo::PlaneGeo collectionPlane = TPC.LastPlane();
      double planeOrigin[3] = {0.};
      double planeCenter[3] = {0.};
      collectionPlane.LocalToWorld(planeOrigin, planeCenter);
      if (TPC.DriftDistance() > 25)
        abs_X_collection = planeCenter[0];
      _activeBounds[0] = -TMath::Abs(abs_X_collection);
      _activeBounds[1] = TMath::Abs(abs_X_collection);
    } // for all TPC
    // set the guard to true
    _isActiveBoundsInitialised = true;
    return;
  }

  // Get active volume
  double *GeometryHelper::GetActiveVolumeBounds() {
    if (!_isActiveBoundsInitialised)
      std::cerr << "Returning active bound vector NOT initialised." << std::endl;
    return _activeBounds;
  }

  // Set fiducial bounds offset from active volume bounds
  void GeometryHelper::SetFiducialBoundOffset(const double &offset) {
    _fiducialBoundOffset = offset;
    _isFiducialBoundOffsetSet = true;
    return;
  }

  // Initialise fiducial volume
  void GeometryHelper::InitFiducialVolumeBounds() {
    if (!_isActiveBoundsInitialised)
      std::cerr << "Returning WRONG fiducial volume." << std::endl;
    if (!_isFiducialBoundOffsetSet)
      std::cerr << "Returning WRONG fiducial volume." << std::endl;
    for (int i=0;i<6;i++) {
      if (i%2==0) _fiducialBounds[i]=_activeBounds[i]+_fiducialBoundOffset;
      else _fiducialBounds[i]=_activeBounds[i]-_fiducialBoundOffset;
    }
    _isFiducialBoundOffsetSet = true;
    return;
  }

  // Get fiducial volume
  double *GeometryHelper::GetFiducialVolumeBounds() {
    if (!_isFiducialBoundOffsetSet)
      std::cerr << "Returning fiducial bound vector NOT initialised." << std::endl;
    return _fiducialBounds;
  }

  // Check if a point is contained in a general volume
  bool GeometryHelper::IsPointInVolume(double *v, TVector3 const &Point) const {
    return (Point.X() >= v[0] && Point.X() <= v[1]
            && Point.Y() >= v[2] && Point.Y() <= v[3]
            && Point.Z() >= v[4] && Point.Z() <= v[5]);
  }

  // Check if a point is contained in a general volume
  bool GeometryHelper::IsPointInVolume(double *v, double *Point) const {
    const TVector3 point_V(Point[0],Point[1],Point[2]);
    return IsPointInVolume(v, point_V);
  }

  // Set the thickness for the slice around the active Volume
  void GeometryHelper::SetThicknessStartVolume(const double &thickness) {
    _thicknessStartVolume = thickness;
    _isThicknessSet = true;
    return;
  }

  // Check if a point is contained in a slice from the active volume
  bool GeometryHelper::IsPointInSlice(TVector3 const &Point) const {
    if (!_isThicknessSet)
      std::cerr << "Thickness is not set." << std::endl;
    if (!_isActiveBoundsInitialised)
      std::cerr << "Active Bounds not initialised." << std::endl;
    return (   (Point.Y()>=(_activeBounds[3]-thicknessStartVolume) && Point.Y()<=v[3])
            || (Point.X()>=_activeBounds[0] && Point.X()<=(_activeBounds[0]+thicknessStartVolume))
            || (Point.X()<=_activeBounds[1] && Point.X()>=(_activeBounds[1]-thicknessStartVolume))
            || (Point.Z()>=_activeBounds[4] && Point.Z()<=(_activeBounds[4]+thicknessStartVolume))
            || (Point.Z()<=_activeBounds[5] && Point.Z()>=(_activeBounds[5]-thicknessStartVolume)));
  }

  // Check if a point is contained in a slice from the active volume
  bool GeometryHelper::IsPointInSlice(double *Point) const {
    const TVector3 point_V(Point[0],Point[1],Point[2]);
    return IsPointInSlice(point_V);
  }

  // Get the APA boundaries (simple version)
  double *GetAPABoundaries() {
    if (!_isActiveBoundsInitialised)
      std::cerr << "Active Bounds not initialised." << std::endl;
    _APABoundaries[0] = _activeBounds[5]/3.;
    _APABoundaries[1] = _activeBounds[5]*2/3.;
  }

} // end of namespace stoppingcosmicmuonselection

#endif
