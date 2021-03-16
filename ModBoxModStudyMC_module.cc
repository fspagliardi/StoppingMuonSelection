///////////////////////////////////////////////////////////////////////
// Class:       ModBoxModStudyMC
// Plugin Type: ******
// File:        ModBoxModStudyMC_module.cc
////////////////////////////////////////////////////////////////////////

#include "ModBoxModStudyMC.h"

namespace stoppingcosmicmuonselection {

  void ModBoxModStudyMC::analyze(art::Event const &evt)
  {
    // increase counter and store event number
    fEvNumber = evt.id().event();
    std::cout << "ModBoxModStudyMC_module is on event: " << fEvNumber << std::endl;
    mf::LogVerbatim("ModBoxModStudyMC") << "ModBoxModStudyMC module on event " << fEvNumber;

    // Set the calibration helper.
    calibHelper.Set(evt);

    // Get handles
    // trackHandle is art::ValidHandle<std::vector<recob::Track>>
    auto const pfparticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
    auto const &recoParticles = *pfparticleHandle;
    auto const spacePointHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointTag);
    const std::vector<recob::SpacePoint> spacePoints = *spacePointHandle;

    auto const hitHandle = evt.getValidHandle<std::vector<recob::Hit>>("hitpdune");
    auto const &allHits = *hitHandle;

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
      fHitAmpl.clear();
      fHitRMS.clear();
      fEfX.clear();
      fEfY.clear();
      fEfZ.clear();
      fEfield.clear();


      // Get the PFParticle
      const recob::PFParticle &thisParticle = recoParticles[p];

      // Only consider primary particles
      if (!thisParticle.IsPrimary()) continue;

      // Skip if the PFParticle is not track-like
      if (!selectorAlg.IsPFParticleATrack(evt,thisParticle)) continue;
      counter_total_number_tracks++;

      //
      //      Make Selection
      //
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

      // Look for and skip track with Michel attached.
      // Get the CNN tagging results.
      size_t trackIndex = hitHelper.GetTrackIndex(track,tracklist);
      auto const &trackHits = hitHelper.GetArtPtrToHitVect(fmht,trackIndex);
      const artPtrHitVec &hitsOnCollection = hitHelper.GetHitsOnAPlane(2,trackHits);
      if (hitsOnCollection.size()==0) continue;

      const size_t &hitIndex = hitHelper.GetIndexClosestHitToPoint(selectorAlg.GetTrackProperties().recoStartPoint,hitsOnCollection,fmthm,tracklist,trackIndex);
      HitPlaneAlg hitPlaneAlg(trackHits,hitIndex,2,selectorAlg.GetTrackProperties().trackT0,clockData,detProp);
      anab::MVAReader<recob::Hit,4> hitResults(evt, fNNetTag);
      if (hitPlaneAlg.AreThereMichelHits(hitResults,0.7,0.5)) continue;

      // Let's go to the Calorimetry. Need to set it for this track first.
      caloHelper.Set(thisParticle,evt,2);
      // Fill the histos
      caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR);
      caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR_TP075,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);
      fdQdx = caloHelper.GetdQdx();
      fDriftTime = caloHelper.GetDriftTime();
      fLifeTimeCorr = caloHelper.GetCorrFactor();
      fResRange = caloHelper.GetResRangeOrdered();
      fTrackPitch = caloHelper.GetTrackPitch();
      fHitX = caloHelper.GetHitX();
      fHitY = caloHelper.GetHitY();
      fHitZ = caloHelper.GetHitZ();
      fPhis = calibHelper.PitchFieldAngle(fHitX, fHitY, fHitZ);
      std::vector<size_t> hitIndeces = caloHelper.GetHitIndex();
      //double xxx = detprop->ConvertTicksToX(allHits[hitIndeces[4]].PeakTime(),allHits[hitIndeces[4]].WireID().Plane, allHits[hitIndeces[4]].WireID().TPC, allHits[hitIndeces[4]].WireID().Cryostat);
      //std::cout << "X: " << fHitX[4] << " Time: " << allHits[hitIndeces[4]].PeakTime() << " Converted: " << xxx << std::endl;
      for (size_t i=0; i<hitIndeces.size();i++) {
        double hitAmpl = allHits[hitIndeces[i]].PeakAmplitude();
        double hitRMS = allHits[hitIndeces[i]].RMS();
        fHitAmpl.push_back(hitAmpl);
        fHitRMS.push_back(hitRMS);
      }

      // Fill the bit for the electric field.
      std::vector<TVector3> electric_field = calibHelper.GetHitPosField(fHitX, fHitY, fHitZ);
      for (size_t i=0; i<electric_field.size(); i++) {
        fEfX.push_back(electric_field.at(i).X());
        fEfY.push_back(electric_field.at(i).Y());
        fEfZ.push_back(electric_field.at(i).Z());
        fEfield.push_back(electric_field.at(i).Mag());
      }

      // Correction for T0.
      // Get rid of cathode crossers for the time being.
      // if ((selectorAlg.GetTrackProperties().trueStartPoint.X()*selectorAlg.GetTrackProperties().trueEndPoint.X()<0) || selectorAlg.GetTrackProperties().trackT0 != INV_DBL)
      //   continue;
      // double x_shift = selectorAlg.GetTrackProperties().trueStartPoint.X() - selectorAlg.GetTrackProperties().recoStartPoint.X();
      // // Correct hit position.
      // for (size_t i=0; i < fHitX.size(); i++) {
      //   fHitX[i] = fHitX[i] + x_shift;
      // }


      // Get Calibration correction factors.
      // TODO: If track is anode crosser, correct x position
      fYZcalibFactor = calibHelper.GetYZCorr_V(fHitX, fHitY, fHitZ);
      fXcalibFactor = calibHelper.GetXCorr_V(fHitX);

      // Correct start and end point.
      sceHelper = new SceHelper(detProp);
      TVector3 recoStartPoint_corr = sceHelper->GetCorrectedPos(TVector3(fStartX, fStartY, fStartZ));
      TVector3 recoEndPoint_corr = sceHelper->GetCorrectedPos(TVector3(fEndX, fEndY, fEndZ));

      fEndX_corr = recoEndPoint_corr.X();
      fEndY_corr = recoEndPoint_corr.Y();
      fEndZ_corr = recoEndPoint_corr.Z();
      fStartX_corr = recoStartPoint_corr.X();
      fStartY_corr = recoStartPoint_corr.Y();
      fStartZ_corr = recoStartPoint_corr.Z();

      // Fill TTree
      std::cout << "*** Adding track..." << std::endl;
      fTrackTree->Fill();

      // Get rid of tracks with weird stuff.
      if (fEndX_corr<-500 or fEndY_corr<-10 or fEndZ_corr<-10 or fStartX_corr<-500 or fStartY_corr<-10 or fStartZ_corr<-10) continue;
      for (size_t i=0; i<fdQdx.size(); i++) {
        if (fResRange.at(i) > 60) {
          h_hitYZ->Fill(fHitZ[i], fHitY[i]);
          h_hitXZ->Fill(fHitZ[i], fHitX[i]);
          h_hitXY->Fill(fHitX[i], fHitY[i]);
        }
      }

    } // end of loop over PFParticles

  } // end of analyzer

} // namespace

DEFINE_ART_MODULE(stoppingcosmicmuonselection::ModBoxModStudyMC)
