/***
  Class containing algorithms for stopping muons selection.

*/
#ifndef CUT_CHECK_HELPER_CXX
#define CUT_CHECK_HELPER_CXX

#include "StoppingMuonSelection/CutCheck/CutCheckHelper.h"

namespace stoppingcosmicmuonselection {

  CutCheckHelper::CutCheckHelper() {
  }

  CutCheckHelper::~CutCheckHelper() {
  }

  // Apply cuts ignoring the specified one.
  void CutCheckHelper::ApplyCutsCathode(TH1 *histo, TH1 *histo_signal, const std::string &excludeCut,
                                        art::Event const &evt, const std::vector<recob::PFParticle> &particles,
                                        const std::vector<recob::SpacePoint> &spacePoints) {

    std::cout << "\tApplying cuts excluding " << excludeCut << std::endl;

    for (unsigned int p = 0; p < particles.size(); p++) {

      // Prepare the selector to digest a new PFParticle
      selectorAlg.Reset();

      // Get the PFParticle
      const recob::PFParticle &thisParticle = particles[p];

      // Only consider primary particles
      if (!thisParticle.IsPrimary()) continue;

      // Skip if the PFParticle is not track-like
      if (!selectorAlg.IsPFParticleATrack(evt,thisParticle)) continue;

      // If this is MC we want that the PFParticle is matched to a cosmic MCParticle
      if (!evt.isRealData() && !selectorAlg.IsTrackMatchedToTrueCosmicTrack(evt,thisParticle))
        continue;

      // Run the selection.
      if (!selectorAlg.NMinus1Cathode(excludeCut, evt, thisParticle)) continue;

      // Check space points.
      const recob::Track &track = selectorAlg.GetTrackFromPFParticle(evt,thisParticle);
      if (!spAlg.IsGoodTrack(track,spacePoints,selectorAlg.GetTrackProperties())) continue;

      const TVector3 &recoStartPoint = selectorAlg.GetTrackProperties().recoStartPoint;
      const TVector3 &recoEndPoint = selectorAlg.GetTrackProperties().recoEndPoint;
      const double &minHitPeakTime = selectorAlg.GetTrackProperties().minHitPeakTime;
      const double &maxHitPeakTime = selectorAlg.GetTrackProperties().maxHitPeakTime;

      if (excludeCut=="thicknessStartVolume") {
        std::string title = histo->GetTitle();
        if (title.find("X")!=std::string::npos)
          histo->Fill(recoStartPoint.X());
        else if (title.find("Y")!=std::string::npos)
          histo->Fill(recoStartPoint.Y());
        else if (title.find("Z")!=std::string::npos)
          histo->Fill(recoStartPoint.Z());
      }
      else if (excludeCut=="cutMinHitPeakTime")
        histo->Fill(minHitPeakTime);
      else if (excludeCut=="cutMaxHitPeakTime")
        histo->Fill(maxHitPeakTime);
      else if (excludeCut == "distanceFiducialVolumeX")
        histo->Fill(recoEndPoint.X());
      else if (excludeCut == "distanceFiducialVolumeY")
        histo->Fill(recoEndPoint.Y());
      else if (excludeCut == "distanceFiducialVolumeZ")
        histo->Fill(recoEndPoint.Z());

      if (!evt.isRealData() && selectorAlg.IsTrueParticleACathodeCrossingStoppingMuon(evt, thisParticle)) {
        if (excludeCut=="thicknessStartVolume") {
          std::string title = histo_signal->GetTitle();
          if (title.find("X")!=std::string::npos)
            histo_signal->Fill(recoStartPoint.X());
          else if (title.find("Y")!=std::string::npos)
            histo_signal->Fill(recoStartPoint.Y());
          else if (title.find("Z")!=std::string::npos)
            histo_signal->Fill(recoStartPoint.Z());
        }
        else if (excludeCut=="cutMinHitPeakTime")
          histo_signal->Fill(minHitPeakTime);
        else if (excludeCut=="cutMaxHitPeakTime")
          histo_signal->Fill(maxHitPeakTime);
        else if (excludeCut == "distanceFiducialVolumeX")
          histo_signal->Fill(recoEndPoint.X());
        else if (excludeCut == "distanceFiducialVolumeY")
          histo_signal->Fill(recoEndPoint.Y());
        else if (excludeCut == "distanceFiducialVolumeZ")
          histo_signal->Fill(recoEndPoint.Z());
      }

      if (excludeCut == "complete") {
        double _trackPitch = 0.75;
        double _trackPitchTolerance = 0.1;
        caloHelper.Set(thisParticle,evt,2);
        if (caloHelper.IsValid()) {
          caloHelper.FillHisto_dQdxVsRR((TH2D *)histo);
          caloHelper.FillHisto_dQdxVsRR((TH2D *)histo_signal,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);
        }
      }

    } // end loop on particles.

  }

  // Apply cuts ignoring the specified one.
  void CutCheckHelper::ApplyCutsAnode(TH1 *histo, TH1 *histo_signal,
                                      const std::string &excludeCut, art::Event const &evt,
                                      const std::vector<recob::PFParticle> &particles,
                                      const std::vector<recob::SpacePoint> &spacePoints) {

    std::cout << "\tApplying cuts excluding " << excludeCut << std::endl;

    for (unsigned int p = 0; p < particles.size(); p++) {

      // Prepare the selector to digest a new PFParticle
      selectorAlg.Reset();

      // Get the PFParticle
      const recob::PFParticle &thisParticle = particles[p];

      // Only consider primary particles
      if (!thisParticle.IsPrimary()) continue;

      // Skip if the PFParticle is not track-like
      if (!selectorAlg.IsPFParticleATrack(evt,thisParticle)) continue;

      // If this is MC we want that the PFParticle is matched to a cosmic MCParticle
      if (!evt.isRealData() && !selectorAlg.IsTrackMatchedToTrueCosmicTrack(evt,thisParticle))
        continue;

      // Run the selection.
      if (!selectorAlg.NMinus1Anode(excludeCut, evt, thisParticle)) continue;

      // Check space points.
      const recob::Track &track = selectorAlg.GetTrackFromPFParticle(evt,thisParticle);
      if (!spAlg.IsGoodTrack(track,spacePoints,selectorAlg.GetTrackProperties())) continue;

      const TVector3 &recoStartPoint = selectorAlg.GetTrackProperties().recoStartPoint;
      const TVector3 &recoEndPoint = selectorAlg.GetTrackProperties().recoEndPoint;
      const double &minHitPeakTime = selectorAlg.GetTrackProperties().minHitPeakTime;
      const double &maxHitPeakTime = selectorAlg.GetTrackProperties().maxHitPeakTime;

      if (excludeCut=="offsetYStartPoint")
        histo->Fill(recoStartPoint.Y());
      else if (excludeCut=="offsetZStartPoint")
        histo->Fill(recoStartPoint.Z());
      else if (excludeCut=="cutMinHitPeakTime")
        histo->Fill(minHitPeakTime);
      else if (excludeCut=="cutMaxHitPeakTime")
        histo->Fill(maxHitPeakTime);
      else if (excludeCut == "distanceFiducialVolumeX")
        histo->Fill(recoEndPoint.X());
      else if (excludeCut == "distanceFiducialVolumeY")
        histo->Fill(recoEndPoint.Y());
      else if (excludeCut == "distanceFiducialVolumeZ")
        histo->Fill(recoEndPoint.Z());

      if (!evt.isRealData() && selectorAlg.IsTrueParticleAnAnodeCrossingStoppingMuon(evt, thisParticle)) {
        if (excludeCut=="offsetYStartPoint")
          histo->Fill(recoStartPoint.Y());
        else if (excludeCut=="offsetZStartPoint")
          histo->Fill(recoStartPoint.Z());
        else if (excludeCut=="cutMinHitPeakTime")
          histo_signal->Fill(minHitPeakTime);
        else if (excludeCut=="cutMaxHitPeakTime")
          histo_signal->Fill(maxHitPeakTime);
        else if (excludeCut == "distanceFiducialVolumeX")
          histo_signal->Fill(recoEndPoint.X());
        else if (excludeCut == "distanceFiducialVolumeY")
          histo_signal->Fill(recoEndPoint.Y());
        else if (excludeCut == "distanceFiducialVolumeZ")
          histo_signal->Fill(recoEndPoint.Z());
      }

      if (excludeCut == "complete") {
        double _trackPitch = 0.75;
        double _trackPitchTolerance = 0.1;
        caloHelper.Set(thisParticle,evt,2);
        caloHelper.FillHisto_dQdxVsRR((TH2D *)histo);
        caloHelper.FillHisto_dQdxVsRR((TH2D *)histo_signal,_trackPitch-_trackPitchTolerance,_trackPitch+_trackPitchTolerance);
      }
    }

  }

  // Fill distribution for every track and for true cathode crossing tracks.
  void CutCheckHelper::FillTruthDistributionCathode(art::Event const &evt, const std::vector<recob::PFParticle> &particles,
                                    TH1D *h_startXPriori, TH1D *h_startX_signalPriori,
                                    TH1D *h_startYPriori, TH1D *h_startY_signalPriori,
                                    TH1D *h_startZPriori, TH1D *h_startZ_signalPriori,
                                    TH1D *h_endXPriori, TH1D *h_endX_signalPriori,
                                    TH1D *h_endYPriori, TH1D *h_endY_signalPriori,
                                    TH1D *h_endZPriori, TH1D *h_endZ_signalPriori,
                                    TH1D *h_minHitPeakTimePriori, TH1D *h_minHitPeakTime_signalPriori,
                                    TH1D *h_maxHitPeakTimePriori, TH1D *h_maxHitPeakTime_signalPriori) {

    std::cout << "\tFill distributions for true cathode crossing muons..." << std::endl;

    for (unsigned int p = 0; p < particles.size(); p++) {

      // Prepare the selector to digest a new PFParticle
      selectorAlg.Reset();

      // Get the PFParticle
      const recob::PFParticle &thisParticle = particles[p];

      // Only consider primary particles
      if (!thisParticle.IsPrimary()) continue;

      // Skip if the PFParticle is not track-like
      if (!selectorAlg.IsPFParticleATrack(evt,thisParticle)) continue;

      // If this is MC we want that the PFParticle is matched to a cosmic MCParticle
      if (!evt.isRealData() && !selectorAlg.IsTrackMatchedToTrueCosmicTrack(evt,thisParticle))
        continue;

      // Run the selection but just to fill the track properties.
      selectorAlg.IsStoppingCathodeCrosser(evt, thisParticle);
      const TVector3 &recoStartPoint = selectorAlg.GetTrackProperties().recoStartPoint;
      const TVector3 &recoEndPoint = selectorAlg.GetTrackProperties().recoEndPoint;
      const double &minHitPeakTime = selectorAlg.GetTrackProperties().minHitPeakTime;
      const double &maxHitPeakTime = selectorAlg.GetTrackProperties().maxHitPeakTime;

      if (!evt.isRealData() && selectorAlg.IsTrueParticleACathodeCrossingStoppingMuon(evt, thisParticle)) {
        h_startX_signalPriori->Fill(recoStartPoint.X());
        h_startY_signalPriori->Fill(recoStartPoint.Y());
        h_startZ_signalPriori->Fill(recoStartPoint.Z());
        h_endX_signalPriori->Fill(recoEndPoint.X());
        h_endY_signalPriori->Fill(recoEndPoint.Y());
        h_endZ_signalPriori->Fill(recoEndPoint.Z());
        h_minHitPeakTime_signalPriori->Fill(minHitPeakTime);
        h_maxHitPeakTime_signalPriori->Fill(maxHitPeakTime);
      }
      h_startXPriori->Fill(recoStartPoint.X());
      h_startYPriori->Fill(recoStartPoint.Y());
      h_startZPriori->Fill(recoStartPoint.Z());
      h_endXPriori->Fill(recoEndPoint.X());
      h_endYPriori->Fill(recoEndPoint.Y());
      h_endZPriori->Fill(recoEndPoint.Z());
      h_minHitPeakTimePriori->Fill(minHitPeakTime);
      h_maxHitPeakTimePriori->Fill(maxHitPeakTime);

    } // loop on particles.

  }

  // Fill distribution for every track and for true cathode crossing tracks.
  void CutCheckHelper::FillTruthDistributionAnode(art::Event const &evt, const std::vector<recob::PFParticle> &particles,
                                  TH1D *h_startYPriori, TH1D *h_startY_signalPriori,
                                  TH1D *h_startZPriori, TH1D *h_startZ_signalPriori,
                                  TH1D *h_endYPriori, TH1D *h_endY_signalPriori,
                                  TH1D *h_endZPriori, TH1D *h_endZ_signalPriori,
                                  TH1D *h_minHitPeakTimePriori, TH1D *h_minHitPeakTime_signalPriori,
                                  TH1D *h_maxHitPeakTimePriori, TH1D *h_maxHitPeakTime_signalPriori) {

    std::cout << "\tFill distributions for true anode crossing muons..." << std::endl;

    for (unsigned int p = 0; p < particles.size(); p++) {

      // Prepare the selector to digest a new PFParticle
      selectorAlg.Reset();

      // Get the PFParticle
      const recob::PFParticle &thisParticle = particles[p];

      // Only consider primary particles
      if (!thisParticle.IsPrimary()) continue;

      // Skip if the PFParticle is not track-like
      if (!selectorAlg.IsPFParticleATrack(evt,thisParticle)) continue;

      // If this is MC we want that the PFParticle is matched to a cosmic MCParticle
      if (!evt.isRealData() && !selectorAlg.IsTrackMatchedToTrueCosmicTrack(evt,thisParticle))
      continue;

      // Run the selection but just to fill the track properties.
      selectorAlg.IsStoppingAnodeCrosser(evt, thisParticle);
      const TVector3 &recoStartPoint = selectorAlg.GetTrackProperties().recoStartPoint;
      const TVector3 &recoEndPoint = selectorAlg.GetTrackProperties().recoEndPoint;
      const double &minHitPeakTime = selectorAlg.GetTrackProperties().minHitPeakTime;
      const double &maxHitPeakTime = selectorAlg.GetTrackProperties().maxHitPeakTime;

      if (!evt.isRealData() && selectorAlg.IsTrueParticleAnAnodeCrossingStoppingMuon(evt, thisParticle)) {
        h_startY_signalPriori->Fill(recoStartPoint.Y());
        h_startZ_signalPriori->Fill(recoStartPoint.Z());
        h_endY_signalPriori->Fill(recoEndPoint.Y());
        h_endZ_signalPriori->Fill(recoEndPoint.Z());
        h_minHitPeakTime_signalPriori->Fill(minHitPeakTime);
        h_maxHitPeakTime_signalPriori->Fill(maxHitPeakTime);
      }
      h_startYPriori->Fill(recoStartPoint.Y());
      h_startZPriori->Fill(recoStartPoint.Z());
      h_endYPriori->Fill(recoEndPoint.Y());
      h_endZPriori->Fill(recoEndPoint.Z());
      h_minHitPeakTimePriori->Fill(minHitPeakTime);
      h_maxHitPeakTimePriori->Fill(maxHitPeakTime);

    }
  }

  // Configure the selector.
  void CutCheckHelper::reconfigure(fhicl::ParameterSet const &p) {

    selectorAlg.reconfigure(p.get<fhicl::ParameterSet>("StoppingMuonSelectionAlg"));
    caloHelper.reconfigure(p.get<fhicl::ParameterSet>("CalorimetryHelper"));
    spAlg.reconfigure(p.get<fhicl::ParameterSet>("SpacePointAlg"));

  }



} // end of namespace stoppingcosmicmuonselection

#endif
