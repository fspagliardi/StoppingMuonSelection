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

    // Configure the selector.
    void reconfigure(fhicl::ParameterSet const &p);

  private:

    StoppingMuonSelectionAlg selectorAlg; // need configuration
    CalorimetryHelper        caloHelper;   // need configuration
  };
}

#endif
