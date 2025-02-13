///////////////////////////////////////////////////////////////////////
// Class:       SelectionStudyProd4
// Plugin Type: ******
// File:        SelectionStudyProd4_module.cc
////////////////////////////////////////////////////////////////////////

#include "SelectionStudyProd4.h"

namespace stoppingcosmicmuonselection {

  void SelectionStudyProd4::analyze(art::Event const &evt)
  {
    // increase counter and store event number
    counter_total_number_events++;
    fEvNumber = evt.id().event();
    std::cout << "SelectionStudyProd4_module is on event: " << fEvNumber <<  " run: " << evt.id().run() << std::endl;
    mf::LogVerbatim("SelectionStudyProd4") << "SelectionStudyProd4 module on event " << fEvNumber;

    // Get handles
    // trackHandle is art::ValidHandle<std::vector<recob::Track>>
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle; // to use with getByLabel to check it's valid
    //auto const pfparticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
    evt.getByLabel(fPFParticleTag, pfparticleHandle);
    if (!pfparticleHandle.isValid()) return;
    auto const &recoParticles = *pfparticleHandle;
    auto const spacePointHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointTag);
    const std::vector<recob::SpacePoint> spacePoints = *spacePointHandle;
    
    // Add handle for clock and detector properties.
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    
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

      fIsRecoSelectedCathodeCrosser = false;
      fIsRecoSelectedAnodeCrosser = false;
      fIsTrueSelectedCathodeCrosser = false;
      fIsTrueSelectedAnodeCrosser = false;

      // Prepare the selector to digest a new PFParticle
      selectorAlg.Reset();

      // Get the PFParticle
      const recob::PFParticle &thisParticle = recoParticles[p];

      // Only consider primary particles
      if (!thisParticle.IsPrimary()) continue;

      // Skip if the PFParticle is not track-like
      if (!selectorAlg.IsPFParticleATrack(evt,thisParticle)) continue;
      counter_total_number_tracks++;

      // If this is MC we want that the PFParticle is matched to a cosmic MCParticle
      if (!evt.isRealData() && !selectorAlg.IsTrackMatchedToTrueCosmicTrack(evt,thisParticle))
        continue;

      // Check if this PFParticle is a stopping muon.
      //      !!!SELECTION STEP!!!
      //
      if (_selectCC && selectorAlg.IsStoppingCathodeCrosser(evt,thisParticle))
        fIsRecoSelectedCathodeCrosser = true;
      else if (_selectAC && selectorAlg.IsStoppingAnodeCrosser(evt,thisParticle))
        fIsRecoSelectedAnodeCrosser = true;
      else
        continue;

      // Check if the track is missing some space points (need to get
      // an handle on the track)
      const recob::Track &track = selectorAlg.GetTrackFromPFParticle(evt,thisParticle);
      if(!spAlg.IsGoodTrack(track,spacePoints,selectorAlg.GetTrackProperties())) {
        std::cout << "Space point alg: " << "TrackID: " << track.ID() << " not accepted." << std::endl;
        continue;
      }

      // Check if the matched PFParticle is a true stopping muon
      if (!evt.isRealData() && fIsRecoSelectedCathodeCrosser)
        fIsTrueSelectedCathodeCrosser = selectorAlg.IsTrueParticleACathodeCrossingStoppingMuon(evt,thisParticle);
      else if (!evt.isRealData() && fIsRecoSelectedAnodeCrosser)
        fIsTrueSelectedAnodeCrosser = selectorAlg.IsTrueParticleAnAnodeCrossingStoppingMuon(evt,thisParticle);

      std::cout << "**************************" << std::endl;
      std::cout << "Track accepted." << std::endl;
      if (fIsRecoSelectedCathodeCrosser)
        std::cout << "Track is a CATHODE crosser." << std::endl;
      else if (fIsRecoSelectedAnodeCrosser) {
        std::cout << "Track is an ANODE crosser. Selected by Pandora? " << selectorAlg.GetTrackProperties().isAnodeCrosserPandora << std::endl;
      }
      std::cout << "Event: " << selectorAlg.GetTrackProperties().evNumber << std::endl;
      std::cout << "trackID: " << selectorAlg.GetTrackProperties().trackID << std::endl;

      // Updating variables to be stored in TTree
      UpdateTTreeVariableWithTrackProperties(selectorAlg.GetTrackProperties());

      // Let's go to the Calorimetry. Need to set it for this track first.
      caloHelper.Set(thisParticle,evt,2);
      // Fill the histos
      caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR);
      caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR_TP075,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);
      //caloHelper.FillHisto_dQdxVsRR_LTCorr(h_dQdxVsRR_LTCorr);
      //caloHelper.FillHisto_dQdxVsRR_LTCorr(h_dQdxVsRR_TP075_LTCorr,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);
      //caloHelper.FillHisto_dQdEVsRR_LTCorr_MC(h_dQdEVsRR_TP075_LTCorr_MC,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);
      //caloHelper.FillHisto_dQdEVsRR_LTCorr_LV(h_dQdEVsRR_TP075_LTCorr_LV,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);

      size_t trackIndex = hitHelper.GetTrackIndex(track,tracklist);
      auto const &trackHits = hitHelper.GetArtPtrToHitVect(fmht,trackIndex);
      size_t numbMichelLikeHits = 0;

      // Init HitPlaneAlg.
      const artPtrHitVec &hitsOnCollection = hitHelper.GetHitsOnAPlane(2,trackHits);
      std::cout << "Hits on collection size: " << hitsOnCollection.size() << std::endl;
      const size_t &hitIndex = hitHelper.GetIndexClosestHitToPoint(selectorAlg.GetTrackProperties().recoStartPoint,hitsOnCollection,fmthm,tracklist,trackIndex);
      HitPlaneAlg hitPlaneAlg(trackHits,hitIndex,2,selectorAlg.GetTrackProperties().trackT0,clockData,detProp);
      // Get the vectors.
      const std::vector<double> &WireIDs = hitPlaneAlg.GetOrderedWireNumb();
      const std::vector<double> &Qs = hitPlaneAlg.GetOrderedQ();
      const std::vector<double> &Dqds = hitPlaneAlg.GetOrderedDqds();
      const std::vector<double> &QsSmooth = hitPlaneAlg.Smoother(Qs,_numberNeighbors);
      const std::vector<double> &DqdsSmooth = hitPlaneAlg.Smoother(Dqds,_numberNeighbors);
      const std::vector<double> &LocalLin = hitPlaneAlg.CalculateLocalLinearity(_numberNeighbors);
      std::cout << "Ordered hit size: " << hitPlaneAlg.GetOrderedHitVec().size() << std::endl;

      for (const art::Ptr<recob::Hit> &hitp : trackHits) {
        if (hitp->WireID().Plane != 2) continue;
        if (evt.isRealData()) continue;
        if (hitHelper.IsHitMichelLike(hitp,selectorAlg.GetTrackProperties().recoEndPoint,fmthm,tracklist,trackIndex,clockData))
          numbMichelLikeHits++;
      }

      // Get the CNN tagging results.
      anab::MVAReader<recob::Hit,4> hitResults(evt, fNNetTag);
      // Store vector of ordered scores.
      std::vector<double> scores = cnnHelper.GetScoreVector(hitResults,hitPlaneAlg.GetOrderedHitVec());
      cnnHelper.FillHitScoreGraph2D(fg_imageScore, hitResults, hitPlaneAlg.GetOrderedHitVec());
      if (!evt.isRealData()) {
        const artPtrHitVec &michelLikeHits = hitHelper.GetMichelLikeHits(hitPlaneAlg.GetOrderedHitVec(),selectorAlg.GetTrackProperties().recoEndPoint,fmthm,tracklist,trackIndex,clockData);
        const artPtrHitVec &muonLikeHits = hitHelper.GetMuonLikeHits(hitPlaneAlg.GetOrderedHitVec(),selectorAlg.GetTrackProperties().recoEndPoint,fmthm,tracklist,trackIndex,clockData);
        f_michelHitsMichelScore = cnnHelper.GetScoreVector(hitResults, michelLikeHits);
        f_muonHitsMichelScore = cnnHelper.GetScoreVector(hitResults, muonLikeHits);
      }
      const artPtrHitVec &hitsNoMichel = hitPlaneAlg.GetHitVecNoMichel(hitResults,_michelScoreThreshold,_michelScoreThresholdAvg);

      if (numbMichelLikeHits > _minNumbMichelLikeHit && !evt.isRealData()) {
        hitHelper.FillTrackGraph2D(fg_imageCollection,hitPlaneAlg.GetOrderedHitVec(),
                                   selectorAlg.GetTrackProperties().recoEndPoint,2,selectorAlg.GetTrackProperties().trackT0,
                                   clockData,detProp);
        hitHelper.FillTrackGraph2D(fg_imageCollectionNoMichel,hitsNoMichel,
                                   selectorAlg.GetTrackProperties().recoEndPoint,2,selectorAlg.GetTrackProperties().trackT0,
                                   clockData,detProp);
          
        // std::cout << "Reco End Point: " << selectorAlg.GetTrackProperties().recoEndPoint.X() << " " << selectorAlg.GetTrackProperties().recoEndPoint.Y() << " " << selectorAlg.GetTrackProperties().recoEndPoint.Z() << std::endl;
        // std::cout << "True End Point: " << selectorAlg.GetTrackProperties().trueEndPoint.X() << " " << selectorAlg.GetTrackProperties().trueEndPoint.Y() << " " << selectorAlg.GetTrackProperties().trueEndPoint.Z() << std::endl;
        // std::cout << "Last hit: " << lastHitPos.X() << " " << lastHitPos.Y() << " " << lastHitPos.Z() << std::endl;
        // std::cout << "Wire: " << geoHelper.GetWireNumb(hitsNoMichel2.back()) << " Time: " << hitsNoMichel2.back()->PeakTime() << std::endl;
      }
      else {
        fg_imageCollection->Set(0);
        fg_imageCollectionNoMichel->Set(0);
      }

      // Fill distance end points.
      fDistEndPoint = (selectorAlg.GetTrackProperties().recoEndPoint - selectorAlg.GetTrackProperties().trueEndPoint).Mag();
      if (!hitPlaneAlg.AreThereMichelHits(hitResults,0.7,0.5)) {
        fDistEndPointNoMichel = (selectorAlg.GetTrackProperties().recoEndPoint - selectorAlg.GetTrackProperties().trueEndPoint).Mag();
        caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR_NoMichel);
        caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR_NoMichelTP,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);
      }
      else
        fDistEndPointNoMichel = INV_DBL;

      // Fill the graphs.
      FillTGraph(fg_wireID, WireIDs);
      FillTGraph(fg_Q,Qs);
      FillTGraph(fg_Dqds,Dqds);
      FillTGraph(fg_QSmooth,QsSmooth);
      FillTGraph(fg_DqdsSmooth,DqdsSmooth);
      FillTGraph(fg_LocalLin,LocalLin);
      FillTGraph(fg_CnnScore, scores);

      fh_progressiveDistance->Reset();
      for (const auto &el : hitPlaneAlg.GetDistances()) {
        fh_progressiveDistance->Fill(el);
      }

      // Fill TTree
      fTrackTree->Fill();

    } // end of loop over PFParticles

  } // end of analyzer

} // namespace

DEFINE_ART_MODULE(stoppingcosmicmuonselection::SelectionStudyProd4)
