/***
  Class containing useful functions for geometry.

*/
#ifndef GEOMETRY_HELPER_H
#define GEOMETRY_HELPER_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "TVector3.h"
#include "TMath.h"

namespace stoppingcosmicmuonselection {

  class GeometryHelper {

  public:
    GeometryHelper();
    ~GeometryHelper();

    // Initialise active volume
    void InitActiveVolumeBounds();

    // Get active volume
    double *GetActiveVolumeBounds();

    // Set fiducial bounds offset from active volume bounds
    void SetFiducialBoundOffset(const double &offset);

    // Initialise fiducial volume
    void InitFiducialVolumeBounds();

    // Get fiducial volume
    double *GetFiducialVolumeBounds();

    // Check if a point is contained in a general volume
    bool IsPointInVolume(double *v, TVector3 const &Point) const;
    bool IsPointInVolume(double *v, double *Point) const;

    // Set the thickness for the slice around the active Volume
    void SetThicknessStartVolume(const double &thickness);

    // Check if a point is contained in a slice from the active volume
    bool IsPointInSlice(TVector3 const &Point) const;
    bool IsPointInSlice(double *Point) const;

    // Get the APA boundaries (simple version)
    double *GetAPABoundaries();

  private:
    bool _isActiveBoundsInitialised = false;
    bool _isFiducialBoundsInitialised = false;
    bool _isThicknessSet = false;
    bool _isFiducialBoundOffsetSet = false;
    double _fiducialBoundOffset;
    double _thicknessStartVolume;
    double _activeBounds[6];
    double _fiducialBounds[6];
    double _APABoundaries[2];
  };
}

#endif
