///////////////////////////////////////////////////////////////////////
// Class:       CosmicStudy
// Plugin Type: ******
// File:        CosmicStudy_module.cc
////////////////////////////////////////////////////////////////////////

#include "CosmicStudy.h"

namespace stoppingcosmicmuonselection {

  void CosmicStudy::analyze(art::Event const &evt)
  {
    // increase counter and store event number
    numberOfEvents++;
    std::cout << "CosmicStudy_module is on event: " << evt.id().event() << std::endl;

    // Get handles
    // trackHandle is art::ValidHandle<std::vector<recob::Track>>
    auto const pfparticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>("pandora");
    auto const &recoParticles = *pfparticleHandle;
    // Iterates over the vector of PFParticles
    for (unsigned int p = 0; p < recoParticles.size(); ++p) {

      // Prepare the selector to digest a new PFParticle
      selectorAlg.Reset();

      // Get the PFParticle
      const recob::PFParticle &thisParticle = recoParticles[p];

      // Only consider primary particles
      if (!thisParticle.IsPrimary()) continue;

      // Skip if the PFParticle is not track-like
      if (!selectorAlg.IsPFParticleATrack(evt,thisParticle)) continue;

      // If this is MC we want that the PFParticle is matched to a cosmic MCParticle
      if (!evt.isRealData() && !selectorAlg.IsTrackMatchedToTrueCosmicTrack(evt,thisParticle))
        continue;

      const simb::MCParticle *particleP = truthUtil.GetMCParticleFromPFParticle(thisParticle,evt,"pandora");
      TVector3 trueEndPoint = particleP->EndPosition().Vect();
      if (abs(particleP->PdgCode()) == 13 && geoHelper.IsPointInVolume(geoHelper.GetActiveVolumeBounds(), trueEndPoint))
        numberOfStoppingMuons++;

      // Get T0
      std::vector<anab::T0> T0s = pfpUtil.GetPFParticleT0(thisParticle, evt, "pandora");
      if (T0s.size()==0) continue;
      if (abs(particleP->PdgCode()) == 13)
        numberOfT0TaggedMuons++;

      // Check if the matched PFParticle is a true stopping muon
      if (selectorAlg.IsTrueParticleACathodeCrossingStoppingMuon(evt,thisParticle))
        numberOfStoppingT0CCMuons++;
      else if (selectorAlg.IsTrueParticleAnAnodeCrossingStoppingMuon(evt,thisParticle))
        numberOfStoppingT0ACMuons++;



    } // end of loop over PFParticles

  } // end of analyzer

} // namespace

DEFINE_ART_MODULE(stoppingcosmicmuonselection::CosmicStudy)
