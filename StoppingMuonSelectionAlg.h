/***
  Class containing useful functions for geometry.

*/
#ifndef STOPPING_MUON_SELECTION_ALG_H
#define STOPPING_MUON_SELECTION_ALG_H

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

#include "GeometryHelper.h"
#include "SpacePointAlg.h"

namespace stoppingcosmicmuonselection {

  class StoppingMuonSelectionAlg {

  public:
    StoppingMuonSelectionAlg();
    ~StoppingMuonSelectionAlg();

    //
    bool IsTrackValid();

    // Read parameters from FHICL file
    void reconfigure(fhicl::ParameterSet const &p);


  private:
    bool _isTrackValid;

    // Helpers and algorithms
    GeometryHelper geoHelper;
    SpacePointAlg  spAlg;

  };
}

#endif
