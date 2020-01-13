#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "TVector3.h"

namespace stoppingcosmicmuonselection {

  constexpr int INV_INT = -999;
  constexpr size_t INV_SIZE = 9999999;
  constexpr double INV_DBL = -9999999;

  typedef std::vector<art::Ptr<recob::Hit>> artPtrHitVec;

  struct trackProperties {
    // Reconstructed information
    size_t evNumber;
    double trackT0;
    TVector3 recoStartPoint;
    TVector3 recoEndPoint;
    double theta_xz, theta_yz;
    double minHitPeakTime, maxHitPeakTime;
    double trackLength;
    double trackID;

    // Characterisation
    bool isCathodeCrosser;
    bool isAnodeCrosserPandora;
    bool isAnodeCrosserMine;

    // Truth information
    int pdg;
    TVector3 trueStartPoint;
    TVector3 trueEndPoint;
    double trueStartT, trueEndT;
    double trueTrackID;

    void Reset() {
      isCathodeCrosser = false;
      isAnodeCrosserPandora = false;
      isAnodeCrosserMine = false;
      evNumber = INV_INT;
      trackT0 = INV_DBL;
      recoStartPoint.SetXYZ(INV_DBL,INV_DBL,INV_DBL);
      recoEndPoint.SetXYZ(INV_DBL,INV_DBL,INV_DBL);
      theta_xz = INV_DBL;
      theta_yz = INV_DBL;
      minHitPeakTime = INV_DBL;
      maxHitPeakTime = INV_DBL;
      trackLength = INV_DBL;
      trackID = INV_DBL;
      pdg = INV_INT;
      trueStartPoint.SetXYZ(INV_DBL,INV_DBL,INV_DBL);
      trueEndPoint.SetXYZ(INV_DBL,INV_DBL,INV_DBL);
      trueStartT = INV_DBL;
      trueEndT = INV_DBL;
      trueTrackID = INV_DBL;
    }
  };

}

#endif
