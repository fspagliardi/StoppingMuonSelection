/***
  Class containing useful functions for calibration.

*/
#ifndef CALIBRATION_HELPER_CXX
#define CALIBRATION_HELPER_CXX

#include "CalibrationHelper.h"

namespace stoppingcosmicmuonselection {

  CalibrationHelper::CalibrationHelper() {

  }

  CalibrationHelper::~CalibrationHelper() {

  }

  // Get the histos.
  void CalibrationHelper::Set(art::Event const &evt) {
    std::string filetype;
    if (!evt.isRealData())
      filetype = "sce";
    else {
      size_t runNumber = evt.id().run();
      filetype = "r" + std::to_string(runNumber);
    }
    std::string filenameX = "Xcalo_" + filetype + ".root";
    std::string filenameYZ = "YZcalo_" + filetype + ".root";
    std::cout << "CalibrationHelper.cxx: filenameX = " << filenameX <<std::endl;
    std::cout << "CalibrationHelper.cxx: filenameYZ = " << filenameYZ <<std::endl;
    TFile fileX(filenameX.c_str());
    TFile fileYZ(filenameYZ.c_str());

    h_x = (TH1D*)fileX.Get("dqdx_X_correction_hist_2");
    h_yz_neg = (TH2D*)fileYZ.Get("correction_dqdx_ZvsY_negativeX_hist_2");
    h_yz_pos = (TH2D*)fileYZ.Get("correction_dqdx_ZvsY_positiveX_hist_2");
  }

  // Get X correction factor.
  double CalibrationHelper::GetXCorr(const TVector3 &hitPos) {
    return h_x->GetBinContent(h_x->FindBin(hitPos.X()));
  }

  // Get YZ correction factor.
  double CalibrationHelper::GetYZCorr(const TVector3 &hitPos) {
    double factor = INV_DBL;

    if (hitPos.X() > 0)
      factor = h_yz_pos->GetBinContent(h_yz_pos->FindBin(hitPos.Z(),hitPos.Y()));
    else
      factor = h_yz_neg->GetBinContent(h_yz_neg->FindBin(hitPos.Z(),hitPos.Y()));

    return factor;
  }

  // Get both factors at the same time.
  double CalibrationHelper::GetXYZCorr(const TVector3 &hitPos) {
    return GetXCorr(hitPos)*GetYZCorr(hitPos);
  }

  // Get vector of factors for X.
  std::vector<double> CalibrationHelper::GetXCorr_V(const std::vector<double> &hit_xs) {
    std::vector<double> result;

    for (auto const &x : hit_xs) {
      TVector3 hitPos(x, 0, 0);
      result.push_back(GetXCorr(hitPos));
    }

    return result;
  }

  // Get vector of factors for YZ.
  std::vector<double> CalibrationHelper::GetYZCorr_V(const std::vector<double> &hit_xs, const std::vector<double> &hit_ys, const std::vector<double> &hit_zs) {
    std::vector<double> result;

    for (size_t id = 0; id < hit_xs.size(); id++) {
      TVector3 hitPos(hit_xs[id], hit_ys[id], hit_zs[id]);
      result.push_back(GetYZCorr(hitPos));
    }

    return result;
  }

} // end of namespace stoppingcosmicmuonselection

#endif
