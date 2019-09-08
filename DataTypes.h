#ifndef DATA_TYPES.H
#define DATA_TYPES.H

namespace stoppingcosmicmuonselection {

  double INV_DBL = -9999999;

  struct trackProperties {
    // Reconstructed information
    size_t evNumber;
    double trackT0;
    double recoStartX, recoStartY, recoStartZ;
    double recoEndX, recoEndY, recoEndZ;
    double theta_xz, theta_yz;
    double minHitPeakTime, maxHitPeakTime;
    double trackLength;
    double trackID;

    // Truth information
    int pdg;
    double trueStartX, trueStartY, trueStartZ;
    double trueEndX, trueEndY, trueEndZ;
    double trueStartT, trueEndT;
    double trueTrackID;

    Reset() {
      evNumber = INV_DBL;
      trackT0 = INV_DBL;
      recoStartX = INV_DBL;
      recoStartY = INV_DBL;
      recoStartZ = INV_DBL;
      recoEndX = INV_DBL;
      recoEndY = INV_DBL;
      recoEndZ = INV_DBL;
      theta_xz = INV_DBL;
      theta_yz = INV_DBL;
      minHitPeakTime = INV_DBL;
      maxHitPeakTime = INV_DBL;
      trackLength = INV_DBL;
      trackID = INV_DBL;
      pdg = INV_DBL;
      trueStartX = INV_DBL;
      trueStartY = INV_DBL;
      trueStartZ = INV_DBL;
      trueEndX = INV_DBL;
      trueEndY = INV_DBL;
      trueEndZ = INV_DBL;
      trueStartT = INV_DBL;
      trueEndT = INV_DBL;
      trueTrackID = INV_DBL;
    }
  }

}

#endif
