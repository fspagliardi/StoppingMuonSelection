#ifndef DATA_TYPES.H
#define DATA_TYPES.H

namespace stoppingcosmicmuonselection {

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
  }
}

#endif
