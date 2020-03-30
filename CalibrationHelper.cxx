/***
  Class containing useful functions for calibration.

*/
#ifndef CALIBRATION_HELPER_CXX
#define CALIBRATION_HELPER_CXX

#include "CalibrationHelper.h"

namespace stoppingcosmicmuonselection {

  CalibrationHelper::CalibrationHelper(const std::string &filetype) {
    std::string filenameX = "XCalo" + filetype;
    std::string filenameYZ = "YZCalo" + filetype;
    TFile fileX(filenameX.c_str());
    TFile fileYZ(filenameYZ.c_str())

    h_x = (TH1D*)fileX.Get("dqdx_X_correction_hist_2");
    h_yz_neg = (TH2D*)fileYZ.Get("correction_dqdx_ZvsY_negativeX_hist_2");
    h_yz_pos = (TH2D*)fileYZ.Get("correction_dqdx_ZvsY_positiveX_hist_2");
  }

  CalibrationHelper::~CalibrationHelper() {

  }

} // end of namespace stoppingcosmicmuonselection

#endif
