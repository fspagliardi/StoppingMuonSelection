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

#include "DataTypes.h"

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
    bool IsPointInVolume(double *v, TVector3 const &Point);
    bool IsPointInVolume(double *v, double *Point);

    // Set the thickness for the slice around the active Volume
    void SetThicknessStartVolume(const double &thickness);

    // Check if a point is contained in a slice from the active volume
    bool IsPointInSlice(TVector3 const &Point);
    bool IsPointInSlice(double *Point);

    // Get the APA boundaries (simple version)
    double *GetAPABoundaries();

    // Get the number of wires from one beam side for a given plane
    size_t GetNumberWiresOneSide(const size_t &planeNumber);

    // Constant to add to number of wires.
    size_t GetWireOffset(const unsigned int &hit_tpcid, const size_t &planeNumber);
    size_t GetWireOffset(const art::Ptr<recob::Hit> &hit, const size_t &planeNumber);

    // Arrays with TPC number info
    const unsigned int tpcIndecesBL[3] = {2,6,10};
    const unsigned int tpcIndecesBR[3] = {1,5,9};

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

    // Handle for geometry service
    const geo::GeometryCore *geom = lar::providerFrom<geo::Geometry>();

  };
}

#endif
