///////////////////////////////////////////////////////////////////////
// Class:       ModBoxModStudyAnode
// Plugin Type: ******
// File:        ModBoxModStudyAnode_module.cc
////////////////////////////////////////////////////////////////////////

#include "ModBoxModStudyMCAnode.h"

namespace stoppingcosmicmuonselection {

  void ModBoxModStudyAnode::analyze(art::Event const &evt)
  {
    // increase counter and store event number
    fEvNumber = evt.id().event();
    std::cout << "ModBoxModStudyAnode_module is on event: " << fEvNumber << std::endl;
    mf::LogVerbatim("ModBoxModStudyAnode") << "ModBoxModStudyAnode module on event " << fEvNumber;
    
    if (evt.isRealData()) {
      art::ServiceHandle<calib::LifetimeCalibService> lifetimecalibHandler;
      calib::LifetimeCalibService & lifetimecalibService = *lifetimecalibHandler;
      calib::LifetimeCalib *lifetimecalib = lifetimecalibService.provider();
      fLifetime = lifetimecalib->GetLifetime()*1000.0; // [ms]*1000.0 -> [us]
    }
    else
      fLifetime = 35000; // [ms]*1000.0 -> [us]
    std::cout << "LIFETIME: " << fLifetime << std::endl;
    
    // Timing stuff
    const char * timestamp;
    art::Timestamp ts = evt.time();
    if (ts.timeHigh()==0) {
      TTimeStamp ts2(ts.timeLow());
      timestamp = ts2.AsString();
    }
    else {
      TTimeStamp ts2(ts.timeHigh(), ts.timeLow());
      timestamp = ts2.AsString();
    }
    std::cout << "TIMESTAMP: "  << timestamp << std::endl;

    // Set the calibration helper.
    calibHelper.Set(evt);

    // Get handles
    // trackHandle is art::ValidHandle<std::vector<recob::Track>>
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle; // to use with getByLabel to check it's valid
    //auto const pfparticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
    evt.getByLabel(fPFParticleTag, pfparticleHandle);
    if (!pfparticleHandle.isValid()) return;
    auto const &recoParticles = *pfparticleHandle;
    auto const spacePointHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointTag);
    const std::vector<recob::SpacePoint> spacePoints = *spacePointHandle;

    auto const hitHandle = evt.getValidHandle<std::vector<recob::Hit>>("hitpdune");
    auto const &allHits = *hitHandle;

    // Add handle for clock and detector properties.
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    double driftVelocity = detProp.DriftVelocity()*1e-3;
    std::cout << "Drift velocity: " << driftVelocity << std::endl;

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
      //caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR);
      //caloHelper.FillHisto_dQdxVsRR(h_dQdxVsRR_TP075,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);
      if (fIsRecoSelectedAnodeCrosser && selectorAlg.GetTrackProperties().isAnodeCrosserPandora) {
        fdQdx = caloHelper.GetdQdx();
        fDriftTime = caloHelper.GetDriftTime();
        fResRange = caloHelper.GetResRangeOrdered();
        fTrackPitch = caloHelper.GetTrackPitch();
        fHitX = caloHelper.GetHitX();
        fHitY = caloHelper.GetHitY();
        fHitZ = caloHelper.GetHitZ();
      }
      else if (fIsRecoSelectedAnodeCrosser && selectorAlg.GetTrackProperties().isAnodeCrosserMine) {
        //calibHelper.CorrectXPosition(fHitX,selectorAlg.GetTrackProperties().recoStartPoint.X(),selectorAlg.GetTrackProperties().recoEndPoint.X(),selectorAlg.GetTrackProperties().trackT0);
        const std::vector<std::vector<double>> &myCalo = fixCalo.GetRightCalo(evt,selectorAlg.GetTrackProperties().trackT0,track);
        if (myCalo.size()!=7) {
          std::cout << "Error: The size of the vector myCalo is wrong!" << std::endl;
          continue;
        }
        fdQdx = myCalo.at(0);
        fResRange = myCalo.at(1);
        fTrackPitch = myCalo.at(2);
        fHitX = myCalo.at(3);
        fHitY = myCalo.at(4);
        fHitZ = myCalo.at(5);
        fdEdx = myCalo.at(6);

        // Apply lifetime correction
        for (size_t j=0;j<fdQdx.size();j++) {
          fdQdx[j] = fdQdx[j] * calibHelper.GetLifeTimeCorrFactor(fLifetime, fHitX[j], evt);
          fLifeTimeCorr.push_back(calibHelper.GetLifeTimeCorrFactor(fLifetime, fHitX[j], evt));
          double ltP10 = 0.1*fLifetime + fLifetime;
          double ltM10 = -0.1*fLifetime + fLifetime;
          fLifeTimeCorrP10.push_back(calibHelper.GetLifeTimeCorrFactor(ltP10, fHitX[j], evt));
          fLifeTimeCorrM10.push_back(calibHelper.GetLifeTimeCorrFactor(ltM10, fHitX[j], evt));
        }

        // Order residual range
        std::vector<double> res_vect;
        std::vector<double> fResRange_ord;
        for (size_t i = 0; i < fResRange.size();i++) {
          res_vect.push_back(fResRange[i]);
        }
        fResRange_ord.resize(fResRange.size());
        size_t size = fResRange.size();
        double max = *max_element(res_vect.begin(),res_vect.end());
        for (size_t i = 0; i < fResRange.size();i++) {
          if (fHitY[size-1] < fHitY[0])  {
            if (fResRange[size-1] < fResRange[0])  {
              fResRange_ord[i] = fResRange[i];
            }
            else {
              fResRange_ord[i] = max-fResRange[i];
            }
          }
          else {
            if (fResRange[size-1] < fResRange[0])  {
              fResRange_ord[i] = max-fResRange[i];
            }
            else {
              fResRange_ord[i] = fResRange[i];
            }
          }
        }
        fResRange = fResRange_ord;
      }
      // Fix lifetime
      //caloHelper.LifeTimeCorrNew(fdQdx, fHitX, evt);
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

      // Get Calibration correction factors.
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

DEFINE_ART_MODULE(stoppingcosmicmuonselection::ModBoxModStudyAnode)
