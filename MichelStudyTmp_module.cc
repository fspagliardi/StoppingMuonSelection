///////////////////////////////////////////////////////////////////////
// Class:       MichelStudyTmp
// Plugin Type: ******
// File:        MichelStudyTmp_module.cc
////////////////////////////////////////////////////////////////////////

#include "MichelStudyTmp.h"

namespace stoppingcosmicmuonselection {

  void MichelStudyTmp::analyze(art::Event const &evt)
  {
    // increase counter and store event number
    counter_total_number_events++;
    fEvNumber = evt.id().event();
    std::cout << "MichelStudyTmp_module is on event: " << fEvNumber << std::endl;
    mf::LogVerbatim("MichelStudyTmp") << "MichelStudyTmp module on event " << fEvNumber;

    // Get handles
    // trackHandle is art::ValidHandle<std::vector<recob::Track>>
    auto const pfparticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
    auto const &recoParticles = *pfparticleHandle;
    auto const spacePointHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointTag);
    const std::vector<recob::SpacePoint> spacePoints = *spacePointHandle;
    // Get track handle.
    art::Handle<std::vector<recob::Track>> trackListHandle;
    std::vector<art::Ptr<recob::Track>> tracklist;
    if (evt.getByLabel(fTrackerTag,trackListHandle))
      art::fill_ptr_vector(tracklist, trackListHandle);
    // Get association to metadata and hit.
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackerTag);
    art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackerTag);

    // Iterates over the vector of PFParticles
    for (unsigned int p = 0; p < recoParticles.size(); ++p) {
      // Prepare the selector to digest a new PFParticle
      selectorAlg.Reset();

      // Get the PFParticle
      const recob::PFParticle thisParticle = recoParticles[p];

      // Only consider primary particles
      if (!thisParticle.IsPrimary()) continue;

      // Skip if the PFParticle is not track-like
      if (!selectorAlg.IsPFParticleATrack(evt,thisParticle)) continue;
      counter_total_number_tracks++;

      // If this is MC we want that the PFParticle is matched to a cosmic MCParticle
      if (!evt.isRealData() && !selectorAlg.IsTrackMatchedToTrueCosmicTrack(evt,thisParticle))
        continue;

      // Check if this PFParticle is a stopping cathode crosser
      if (!selectorAlg.IsStoppingCathodeCrosser(evt,thisParticle))
        continue;

      // Check if the track is missing some space points (need to get
      // an handle on the track)
      const recob::Track &track = selectorAlg.GetTrackFromPFParticle(evt,thisParticle);
      spAlg.SetT0(selectorAlg.GetTrackProperties().trackT0);
      if(!spAlg.IsGoodTrack(track,spacePoints))
        continue;

      // At this point we have selected a cathode crossing stopping muon
      fIsRecoSelectedCathodeCrosser = true;

      // Check if the matched PFParticle is a cathode-crossing stopping muon
      fIsTrueSelectedCathodeCrosser = selectorAlg.IsTrueParticleAStoppingMuon(evt,thisParticle);

      std::cout << "**************************" << std::endl;
      std::cout << "Track accepted." << std::endl;
      std::cout << "Event: " << selectorAlg.GetTrackProperties().evNumber << std::endl;
      std::cout << "trackID: " << selectorAlg.GetTrackProperties().trackID << std::endl;

      // Updating variables to be stored in TTree
      UpdateTTreeVariableWithTrackProperties(selectorAlg.GetTrackProperties());

      // Let's go to the Calorimetry. Need to set it for this track first.
      caloHelper.Set(thisParticle,evt);
      // Fill the histos
      caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR,2);
      caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR_TP075,2,0.75-0.1,0.75+0.1);
      caloHelper.FillHisto_dQdxVsRR_LTCorr(h_dQdxVsRR_LTCorr,2);
      caloHelper.FillHisto_dQdxVsRR_LTCorr(h_dQdxVsRR_TP075_LTCorr,2,0.75-0.1,0.75+0.1);

      size_t trackIndex = hitHelper.GetTrackIndex(track,tracklist);
      auto const &trackHits = hitHelper.GetArtPtrToHitVect(fmht,trackIndex);
      size_t numbMichelLikeHits = 0;
      for (const art::Ptr<recob::Hit> &hitp : trackHits) {
        if (hitHelper.IsHitMichelLike(hitp,selectorAlg.GetTrackProperties().recoEndPoint,fmthm,tracklist,trackIndex))
          numbMichelLikeHits++;
      }
      if (numbMichelLikeHits > _minNumbMichelLikeHit)
        hitHelper.FillTrackHitPicture(fh_imageCollection,trackHits,selectorAlg.GetTrackProperties().recoEndPoint,2);
      else
        fh_imageCollection->Reset();

      // Init HitPlaneAlg.
      const artPtrHitVec &hitsOnCollection = hitHelper.GetHitsOnAPlane(2,trackHits);
      const size_t &hitIndex = hitHelper.GetIndexClosestHitToPoint(selectorAlg.GetTrackProperties().recoStartPoint,hitsOnCollection,fmthm,tracklist,trackIndex);
      HitPlaneAlg hitPlaneAlg(trackHits,hitIndex,2);
      // Create the graph
      hitPlaneAlg.FillTGraphDqds(g_Q);
      hitPlaneAlg.FillTGraphDqds(g_Dqds);

      // Fill TTree
      fTrackTree->Fill();

    } // end of loop over PFParticles

  } // end of analyzer

} // namespace

DEFINE_ART_MODULE(stoppingcosmicmuonselection::MichelStudyTmp)
