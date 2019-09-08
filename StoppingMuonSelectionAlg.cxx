/***
  Class containing algorithms for stopping muons selection.

*/
#ifndef STOPPING_MUON_SELECTION_ALG_CXX
#define STOPPING_MUON_SELECTION_ALG_CXX

#include "StoppingMuonSelectionAlg.h"

namespace stoppingcosmicmuonselection {

  StoppingMuonSelectionAlg::StoppingMuonSelectionAlg() {

  }

  StoppingMuonSelectionAlg::~StoppingMuonSelectionAlg() {

  }

  //
  bool IsPFParticleATrack(art::Event const &evt,recob::PFParticle const &thisParticle) {
    const recob::Track *trackP = pfpUtil.GetPFParticleTrack(thisParticle,evt,fPFParticleTag,fTrackerTag);
    if (trackP == nullptr) return false;
    else return true;
  }

  // Read parameters from FHICL file
  void reconfigure(fhicl::ParameterSet const &p) {
    fTrackerTag = p.get<std::string>("TrackerTag");
    fPFParticleTag = p.get<std::string>("PFParticleTag");
    // Set cuts
    // Define cuts
    double length_cutoff_CC = p.get<double>("length_cutoff_CC",100);
    double offsetFiducialBounds_CC = p.get<double>("offsetFiducialBounds_CC",50);
    double thicknessStartVolume_CC = p.get<double>("thicknessStartVolume_CC",40);
    double cutMinHitPeakTime_CC = p.get<double>("cutMinHitPeakTime_CC",500);
    double cutMaxHitPeakTime_CC = p.get<double>("cutMaxHitPeakTime_CC",4800);
    double radiusBrokenTracksSearch_CC = p.get<double>("radiusBrokenTracksSearch_CC",50);
    double cutCosAngleBrokenTracks_CC = p.get<double>("cutCosAngleBrokenTracks_CC",0.995);
    double cutCosAngleAlignment_CC = p.get<double>("cutCosAngleAlignment_CC",0.995);
    double cutContourAPA_CC = p.get<double>("cutContourAPA_CC",10);
  }

  // Determine if the PFParticle is a selected cathode crosser
  bool IsStoppingCathodeCrosser(art::Event const &evt, recob::PFParticle const &thisParticle) {
    // Check that this PFParticle is a track
    if (!IsPFParticleATrack(evt,thisParticle)) return false;
    // In case this event is MC, check if we have a MCParticle from cosmics
    // associated to this PFParticle
    if (!evt.isRealData() && !IsTrackMatchedToTrueCosmicTrack(evt,thisParticle)) return false;
    // Check if it's primary
    if (!thisParticle.IsPrimary()) return false;

    // All good now, proceed with storing info for selection.
    // Get the T0
    std::vector<anab::T0> pfparticleT0s = pfpUtil.GetPFParticleT0(thisParticle,evt,fPFParticleTag);
    if (pfparticleT0s.size() == 0)
      trackT0 = INV_DBL;
    else
      trackT0 = pfparticleT0s[0].Time();
    // Get recob::Track from PFParticle
    const recob::Track &track = GetTrackFromPFParticle(evt,thisParticle);
    recoEndPoint = track.End<TVector3>();
    recoStartPoint.Set(track.LocationAtPoint(track.FirstValidPoint()).X(), track.LocationAtPoint(track.FirstValidPoint()).Y(), track.LocationAtPoint(track.FirstValidPoint()).Z());
    OrderRecoStartEnd(recoStartPoint, recoEndPoint);
    trackLength = track.Length();
    trackID = track.ID();
    // using the ordered start and end points calculate the angles theta_xz and theta_yz
    theta_xz = TMath::RadToDeg() * TMath::ATan2(recoStartPoint.X()-recoEndPoint.X(), recoStartPoint.Z()-recoEndPoint.Z());
    theta_yz = TMath::RadToDeg() * TMath::ATan2(recoStartPoint.Y()-recoEndPoint.Y(), recoStartPoint.Z()-recoEndPoint.Z());
    // Determine min and max hit peak time for this track
    SetMinAndMaxHitPeakTime(evt,thisParticle,minHitPeakTime,maxHitPeakTime);

    // Prepare geometry helper
    geoHelper.SetFiducialBoundOffset(offsetFiducialBounds_CC);
    geoHelper.SetThicknessStartVolume(thicknessStartVolume_CC);
    geoHelper.InitActiveVolumeBounds();

    // Apply cuts with selection with progressive cuts
    if (trackLength < length_cutoff_CC) return false;
    if (TMath::Abs(theta_yz-90)<10 || TMath::Abs(theta_yz+90)<10 || TMath::Abs(theta_xz-90)<10 || TMath::Abs(theta_xz+90)<10) return false;
    bool goodTrack = ((recoStartPoint.X()*recoEndPoint.X()<0)
                     && geoHelper.IsPointInSlice(recoStartPoint)
                     && geoHelper.IsPointInVolume(geoHelper.GetFiducialVolumeBounds(),recoEndPoint));
    if (!goodTrack) return false;
    if (minHitPeakTime <= cutMinHitPeakTime_CC) return false;
    if (maxHitPeakTime >= cutMaxHitPeakTime_CC) return false;


  }

  // For MC events, check if the track is associated to a cosmic track
  bool IsTrackMatchedToTrueCosmicTrack(art::Event const &evt, recob::PFParticle const &thisParticle) {
    const simb::MCParticle *particleP = 0x0;
    if (!evt.isRealData()) {
      particleP = truthUtil.GetMCParticleFromPFParticle(thisParticle,evt,fPFParticleTag);
      if (particleP == 0x0) return false;
      if (pi_serv->TrackIdToMCTruth_P(particleP->TrackId())->Origin() != simb::kCosmicRay))
        return true;
      else
        return false;
    }
    else
      return false;
  }

  // Order reco start and end point based on Y position
  void OrderRecoStartEnd(TVector3 start, TVector3 end) {
    TVector3 prov;
    if (end.Y() > start.Y()) {
      prov = start;
      start = end;
      end = prov;
    }
    return;
  }

  // Get track from PFParticle
  const recob::Track GetTrackFromPFParticle(art::Event const &evt, recob::PFParticle const &thisParticle) {
    const recob::Track *trackP = pfpUtil.GetPFParticleTrack(thisParticle,evt,fPFParticleTag,fTrackerTag);
    return *(trackP);
  }

  // Determine min and max hit peak time for this track
  void SetMinAndMaxHitPeakTime(art::Event const &evt, recob::PFParticle const &thisParticle, duoble &minHitPeakTime, double &maxHitPeakTime) {
    // Get Hits associated with PFParticle
    const std::vector<const recob::Hit*> Hits = pfpUtil.GetPFParticleHits(thisParticle,evt,fPFParticleTag);
    std::vector<double> HitPeakTimes;
    for (unsigned int hitIndex = 0; hitIndex < Hits.size(); ++hitIndex)   {
      HitPeakTimes.push_back(Hits[hitIndex]->PeakTime());
    }
    minHitPeakTime = *(std::min_element(HitPeakTimes.begin(), HitPeakTimes.end()));
    maxHitPeakTime = *(std::max_element(HitPeakTimes.begin(), HitPeakTimes.end()));
  }


} // end of namespace stoppingcosmicmuonselection

#endif
