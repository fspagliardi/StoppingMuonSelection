///////////////////////////////////////////////////////////////////////
// Class:       CutCheck
// Plugin Type: ******
// File:        CutCheck_module.cc
////////////////////////////////////////////////////////////////////////

#include "StoppingMuonSelection/CutCheck/CutCheck.h"

namespace stoppingcosmicmuonselection {

  void CutCheck::analyze(art::Event const &evt)
  {
    evNumber = evt.id().event();
    std::cout << "CutCheck_module is on event: " << evNumber << std::endl;
    mf::LogVerbatim("CutCheck") << "CutCheck module on event " << evNumber;

    // Get handles
    // trackHandle is art::ValidHandle<std::vector<recob::Track>>
    auto const pfparticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
    auto const &recoParticles = *pfparticleHandle;
    auto const spacePointHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointTag);
    const std::vector<recob::SpacePoint> spacePoints = *spacePointHandle;

    if (_selectCC) {
      std::cout << "Analysing cathode-crossers..." << std::endl;
      if (!_runCathodeSimple) {
        cutCheckHelper.ApplyCutsCathode(h_startX, h_startX_signal, "thicknessStartVolume", evt, recoParticles, spacePoints);
    	  cutCheckHelper.ApplyCutsCathode(h_startY, h_startY_signal, "thicknessStartVolume", evt, recoParticles, spacePoints);
    	  cutCheckHelper.ApplyCutsCathode(h_startZ, h_startZ_signal, "thicknessStartVolume", evt, recoParticles, spacePoints);
        cutCheckHelper.ApplyCutsCathode(h_minHitPeakTime, h_minHitPeakTime_signal, "cutMinHitPeakTime", evt, recoParticles, spacePoints);
        cutCheckHelper.ApplyCutsCathode(h_maxHitPeakTime, h_maxHitPeakTime_signal, "cutMaxHitPeakTime", evt, recoParticles, spacePoints);
    	  cutCheckHelper.ApplyCutsCathode(h_endX, h_endX_signal, "distanceFiducialVolumeX", evt, recoParticles, spacePoints);
    	  cutCheckHelper.ApplyCutsCathode(h_endY, h_endY_signal, "distanceFiducialVolumeY", evt, recoParticles, spacePoints);
    	  cutCheckHelper.ApplyCutsCathode(h_endZ, h_endZ_signal, "distanceFiducialVolumeZ", evt, recoParticles, spacePoints);
        cutCheckHelper.ApplyCutsCathode(h_dQdxVsRR, h_dQdxVsRR_TP, "complete", evt, recoParticles, spacePoints);
      }
      else {
        cutCheckHelper.ApplyCutsCathodeSimple(h_startX, h_startX_signal, "thicknessStartVolume", evt, recoParticles);
    	  cutCheckHelper.ApplyCutsCathodeSimple(h_startY, h_startY_signal, "thicknessStartVolume", evt, recoParticles);
    	  cutCheckHelper.ApplyCutsCathodeSimple(h_startZ, h_startZ_signal, "thicknessStartVolume", evt, recoParticles);
        cutCheckHelper.ApplyCutsCathodeSimple(h_minHitPeakTime, h_minHitPeakTime_signal, "cutMinHitPeakTime", evt, recoParticles);
        cutCheckHelper.ApplyCutsCathodeSimple(h_maxHitPeakTime, h_maxHitPeakTime_signal, "cutMaxHitPeakTime", evt, recoParticles);
    	  cutCheckHelper.ApplyCutsCathodeSimple(h_endX, h_endX_signal, "distanceFiducialVolumeX", evt, recoParticles);
    	  cutCheckHelper.ApplyCutsCathodeSimple(h_endY, h_endY_signal, "distanceFiducialVolumeY", evt, recoParticles);
    	  cutCheckHelper.ApplyCutsCathodeSimple(h_endZ, h_endZ_signal, "distanceFiducialVolumeZ", evt, recoParticles);
        cutCheckHelper.ApplyCutsCathodeSimple(h_dQdxVsRR, h_dQdxVsRR_TP, "complete", evt, recoParticles);
      }

      cutCheckHelper.FillTruthDistributionCathode(evt, recoParticles,
                                                  h_startXPriori, h_startX_signalPriori,
                                                  h_startYPriori, h_startY_signalPriori,
                                                  h_startZPriori, h_startZ_signalPriori,
                                                  h_endXPriori, h_endX_signalPriori,
                                                  h_endYPriori, h_endY_signalPriori,
                                                  h_endZPriori, h_endZ_signalPriori,
                                                  h_minHitPeakTimePriori, h_minHitPeakTime_signalPriori,
                                                  h_maxHitPeakTimePriori, h_maxHitPeakTime_signalPriori);

    }
    else if (_selectAC) {
      std::cout << "Analysing anode-crossers..." << std::endl;
      cutCheckHelper.ApplyCutsAnode(h_startY, h_startY_signal, "offsetYStartPoint", evt, recoParticles, spacePoints);
    	cutCheckHelper.ApplyCutsAnode(h_startZ, h_startZ_signal, "offsetZStartPoint", evt, recoParticles, spacePoints);
      cutCheckHelper.ApplyCutsAnode(h_minHitPeakTime, h_minHitPeakTime_signal, "cutMinHitPeakTime", evt, recoParticles, spacePoints);
      cutCheckHelper.ApplyCutsAnode(h_maxHitPeakTime, h_maxHitPeakTime_signal, "cutMaxHitPeakTime", evt, recoParticles, spacePoints);
    	cutCheckHelper.ApplyCutsAnode(h_endX, h_endX_signal, "distanceFiducialVolumeX", evt, recoParticles, spacePoints);
    	cutCheckHelper.ApplyCutsAnode(h_endY, h_endY_signal, "distanceFiducialVolumeY", evt, recoParticles, spacePoints);
    	cutCheckHelper.ApplyCutsAnode(h_endZ, h_endZ_signal, "distanceFiducialVolumeZ", evt, recoParticles, spacePoints);
      cutCheckHelper.ApplyCutsAnode(h_dQdxVsRR, h_dQdxVsRR_TP, "complete", evt, recoParticles, spacePoints);

      cutCheckHelper.FillTruthDistributionAnode(evt, recoParticles,
                                                  h_startYPriori, h_startY_signalPriori,
                                                  h_startZPriori, h_startZ_signalPriori,
                                                  h_endYPriori, h_endY_signalPriori,
                                                  h_endZPriori, h_endZ_signalPriori,
                                                  h_minHitPeakTimePriori, h_minHitPeakTime_signalPriori,
                                                  h_maxHitPeakTimePriori, h_maxHitPeakTime_signalPriori);
    }

    h_events->Fill(0);
  } // end of analyzer

} // namespace

DEFINE_ART_MODULE(stoppingcosmicmuonselection::CutCheck)
