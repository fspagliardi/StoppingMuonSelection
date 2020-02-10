///////////////////////////////////////////////////////////////////////
// Class:       CutCheck
// Plugin Type: ******
// File:        CutCheck_module.cc
////////////////////////////////////////////////////////////////////////

#include "StoppingMuonSelection/CutCheck/CutCheck.h"

namespace stoppingcosmicmuonselection {

  void CutCheck::analyze(art::Event const &evt)
  {
    fEvNumber = evt.id().event();
    std::cout << "CutCheck_module is on event: " << fEvNumber << std::endl;
    mf::LogVerbatim("CutCheck") << "CutCheck module on event " << fEvNumber;

    // Get handles
    // trackHandle is art::ValidHandle<std::vector<recob::Track>>
    auto const pfparticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
    auto const &recoParticles = *pfparticleHandle;

    if (_selectCC) {
      std::cout << "Analysing cathode-crossers..." << std::endl;
      cutCheckHelper.ApplyCutsCathode(h_startX, h_startX_signal, "thicknessStartVolume", evt, recoParticles);
    	cutCheckHelper.ApplyCutsCathode(h_startY, h_startY_signal, "thicknessStartVolume", evt, recoParticles);
    	cutCheckHelper.ApplyCutsCathode(h_startZ, h_startZ_signal, "thicknessStartVolume", evt, recoParticles);
      cutCheckHelper.ApplyCutsCathode(h_minHitPeakTime, h_minHitPeakTime_signal, "cutMinHitPeakTime", evt, recoParticles);
      cutCheckHelper.ApplyCutsCathode(h_maxHitPeakTime, h_maxHitPeakTime_signal, "cutMaxHitPeakTime", evt, recoParticles);
    	cutCheckHelper.ApplyCutsCathode(h_endX, h_endX_signal, "distanceFiducialVolumeX", evt, recoParticles);
    	cutCheckHelper.ApplyCutsCathode(h_endY, h_endY_signal, "distanceFiducialVolumeY", evt, recoParticles);
    	cutCheckHelper.ApplyCutsCathode(h_endZ, h_endZ_signal, "distanceFiducialVolumeZ", evt, recoParticles);
      cutCheckHelper.ApplyCutsCathode(h_dQdxVsRR, h_dQdxVsRR_TP, "complete", evt, recoParticles);
    }


  } // end of analyzer

} // namespace

DEFINE_ART_MODULE(stoppingcosmicmuonselection::CutCheck)
