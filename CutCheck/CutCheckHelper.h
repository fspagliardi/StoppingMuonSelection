/***
  Class containing algorithms for stopping muons selection.

*/
#ifndef CUT_CHECK_HELPER_H
#define CUT_CHECK_HELPER_H

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "larsim/Simulation/LArVoxelData.h"
#include "larsim/Simulation/LArVoxelList.h"
#include "larsim/Simulation/SimListUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"

#include "StoppingMuonSelection/GeometryHelper.h"
#include "StoppingMuonSelection/DataTypes.h"
#include "StoppingMuonSelection/CalorimetryHelper.h"
#include "StoppingMuonSelection/StoppingMuonSelectionAlg.h"

namespace stoppingcosmicmuonselection {

  class CutCheckHelper {

  public:
    CutCheckHelper();
    ~CutCheckHelper();

    // Apply cuts ignoring the specified one.
    void ApplyCutsCathode(TH1 *histo, TH1 *histo_signal, const std::string &excludeCut, art::Event const &evt, const std::vector<recob::PFParticle> &particles);

    // Apply cuts ignoring the specified one.
    void ApplyCutsAnode(TH1 *histo, TH1 *histo_signal, const std::string &excludeCut, art::Event const &evt, const std::vector<recob::PFParticle> &particles);

    // Fill distribution for every track and for true cathode crossing tracks.
    void FillTruthDistributionCathode(art::Event const &evt, const std::vector<recob::PFParticle> &particles,
                                      TH1D *h_startXPriori, TH1D *h_startX_signalPriori,
                                      TH1D *h_startYPriori, TH1D *h_startY_signalPriori,
                                      TH1D *h_startZPriori, TH1D *h_startZ_signalPriori,
                                      TH1D *h_endXPriori, TH1D *h_endX_signalPriori,
                                      TH1D *h_endYPriori, TH1D *h_endY_signalPriori,
                                      TH1D *h_endZPriori, TH1D *h_endZ_signalPriori,
                                      TH1D *h_minHitPeakTimePriori, TH1D *h_minHitPeakTime_signalPriori,
                                      TH1D *h_maxHitPeakTimePriori, TH1D *h_maxHitPeakTime_signalPriori);

    // Configure the selector.
    void reconfigure(fhicl::ParameterSet const &p);

  private:

    StoppingMuonSelectionAlg selectorAlg; // need configuration
    CalorimetryHelper        caloHelper;   // need configuration
  };
}

#endif
