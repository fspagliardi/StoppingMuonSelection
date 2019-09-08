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

  // See if there is track associated to this PFParticle
  bool StoppingMuonSelectionAlg::IsPFParticleATrack(art::Event const &evt,
                                                    recob::PFParticle const &thisParticle) {
    const recob::Track *trackP = pfpUtil.GetPFParticleTrack(thisParticle,evt,fPFParticleTag,fTrackerTag);
    if (trackP == nullptr) return false;
    else return true;
  }

  // Read parameters from FHICL file
  void StoppingMuonSelectionAlg::reconfigure(fhicl::ParameterSet const &p) {
    fTrackerTag = p.get<std::string>("TrackerTag");
    fPFParticleTag = p.get<std::string>("PFParticleTag");
    // Set cuts
    // Define cuts
    length_cutoff_CC = p.get<double>("length_cutoff_CC",100);
    offsetFiducialBounds_CC = p.get<double>("offsetFiducialBounds_CC",50);
    thicknessStartVolume_CC = p.get<double>("thicknessStartVolume_CC",40);
    cutMinHitPeakTime_CC = p.get<double>("cutMinHitPeakTime_CC",500);
    cutMaxHitPeakTime_CC = p.get<double>("cutMaxHitPeakTime_CC",4800);
    radiusBrokenTracksSearch_CC = p.get<double>("radiusBrokenTracksSearch_CC",50);
    cutCosAngleBrokenTracks_CC = p.get<double>("cutCosAngleBrokenTracks_CC",0.995);
    cutCosAngleAlignment_CC = p.get<double>("cutCosAngleAlignment_CC",0.995);
    cutContourAPA_CC = p.get<double>("cutContourAPA_CC",10);
  }

  // Determine if the PFParticle is a selected cathode crosser
  bool StoppingMuonSelectionAlg::IsStoppingCathodeCrosser(art::Event const &evt,
                                                          recob::PFParticle const &thisParticle) {
    // Check that this PFParticle is a track
    if (!IsPFParticleATrack(evt,thisParticle)) return false;
    // In case this event is MC, check if we have a MCParticle from cosmics
    // associated to this PFParticle
    if (!evt.isRealData() && !IsTrackMatchedToTrueCosmicTrack(evt,thisParticle)) return false;
    // Check if it's primary
    if (!thisParticle.IsPrimary()) return false;

    // All good now, proceed with storing info for selection.
    _evNumber = evt.id().event();
    // Get the T0
    std::vector<anab::T0> pfparticleT0s = pfpUtil.GetPFParticleT0(thisParticle,evt,fPFParticleTag);
    if (pfparticleT0s.size() == 0) {
      _trackT0 = INV_DBL;
      return false;
    }
    else
      _trackT0 = pfparticleT0s[0].Time();
    // Get recob::Track from PFParticle
    const recob::Track &track = GetTrackFromPFParticle(evt,thisParticle);
    _recoEndPoint = track.End<TVector3>();
    _recoStartPoint.SetXYZ(track.LocationAtPoint(track.FirstValidPoint()).X(), track.LocationAtPoint(track.FirstValidPoint()).Y(), track.LocationAtPoint(track.FirstValidPoint()).Z());
    OrderRecoStartEnd(_recoStartPoint, _recoEndPoint);
    _trackLength = track.Length();
    _trackID = track.ID();
    // using the ordered start and end points calculate the angles _theta_xz and _theta_yz
    _theta_xz = TMath::RadToDeg() * TMath::ATan2(_recoStartPoint.X()-_recoEndPoint.X(), _recoStartPoint.Z()-_recoEndPoint.Z());
    _theta_yz = TMath::RadToDeg() * TMath::ATan2(_recoStartPoint.Y()-_recoEndPoint.Y(), _recoStartPoint.Z()-_recoEndPoint.Z());
    // Determine min and max hit peak time for this track
    SetMinAndMaxHitPeakTime(evt,thisParticle,_minHitPeakTime,_maxHitPeakTime);

    // Prepare geometry helper
    geoHelper.SetFiducialBoundOffset(offsetFiducialBounds_CC);
    geoHelper.SetThicknessStartVolume(thicknessStartVolume_CC);
    geoHelper.InitActiveVolumeBounds();

    // Apply cuts with selection with progressive cuts
    if (_trackLength < length_cutoff_CC) return false;
    if (TMath::Abs(_theta_yz-90)<10 || TMath::Abs(_theta_yz+90)<10 || TMath::Abs(_theta_xz-90)<10 || TMath::Abs(_theta_xz+90)<10) return false;
    bool goodTrack = ((_recoStartPoint.X()*_recoEndPoint.X()<0)
                     && geoHelper.IsPointInSlice(_recoStartPoint)
                     && geoHelper.IsPointInVolume(geoHelper.GetFiducialVolumeBounds(),_recoEndPoint));
    if (!goodTrack) return false;
    if (_minHitPeakTime <= cutMinHitPeakTime_CC) return false;
    if (_maxHitPeakTime >= cutMaxHitPeakTime_CC) return false;
    if ((TMath::Abs(_recoEndPoint.Z()-geoHelper.GetAPABoundaries()[0])<=cutContourAPA_CC) || (TMath::Abs(_recoEndPoint.Z()-geoHelper.GetAPABoundaries()[1])<=cutContourAPA_CC)) return false;

    // Look for broken tracks
    TVector3 dirFirstTrack = _recoEndPoint - _recoStartPoint;
    bool isBrokenTrack = false;
    const std::vector<recob::PFParticle> & pfparticles = *(evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag));
    for (size_t p=0;p<pfparticles.size();p++) {
      // Get track
      const recob::Track *newTrack = pfpUtil.GetPFParticleTrack(pfparticles[p],evt,fPFParticleTag,fTrackerTag);
      if (newTrack==nullptr) continue;
      if (newTrack->ID()==_trackID) continue;
      size_t fp = newTrack->FirstValidPoint();
      TVector3 recoStartPointSecond(newTrack->LocationAtPoint(fp).X(),newTrack->LocationAtPoint(fp).Y(),newTrack->LocationAtPoint(fp).Z());
      TVector3 recoEndPointSecond = newTrack->End<TVector3>();
      OrderRecoStartEnd(recoStartPointSecond,recoEndPointSecond);
      TVector3 dirSecondTrack = recoEndPointSecond-recoStartPointSecond;
      TVector3 dirHigherTrack, dirLowerTrack, endPointHigherTrack, startPointLowerTrack, endPointLowerTrack;
      if (_recoStartPoint.Y() > recoStartPointSecond.Y()) {
        dirHigherTrack = dirFirstTrack;
        dirLowerTrack = dirSecondTrack;
        startPointLowerTrack = recoStartPointSecond;
        endPointLowerTrack = recoEndPointSecond;
        endPointHigherTrack = _recoEndPoint;
      }
      else
        continue;
      TVector3 middePointLowerTrack = (startPointLowerTrack + endPointLowerTrack) * 0.5;
      TVector3 dirHigherTrack_YZ(0., dirHigherTrack.Y(), dirHigherTrack.Z());
      TVector3 dirJoiningSegment(0., middePointLowerTrack.Y()-endPointHigherTrack.Y(), middePointLowerTrack.Z()-endPointHigherTrack.Z());
      double cosBeta = TMath::Cos(dirHigherTrack_YZ.Angle(dirJoiningSegment));
      double absCosAlpha = TMath::Abs(TMath::Cos(dirFirstTrack.Angle(dirSecondTrack)));

      if ( (absCosAlpha > cutCosAngleBrokenTracks_CC) && (cosBeta >= cutCosAngleAlignment_CC) )
        isBrokenTrack = true;

      double distHigherLower = TMath::Sqrt(TMath::Power(endPointHigherTrack.Y()-startPointLowerTrack.Y(),2) + TMath::Power(endPointHigherTrack.Z()-startPointLowerTrack.Z(),2));
      if (distHigherLower < radiusBrokenTracksSearch_CC) {
        if ( TMath::Abs(TMath::Cos(dirFirstTrack.Angle(dirSecondTrack))) > 0.96)
          isBrokenTrack = true;
      }
      // Stop searching if found one
      if (isBrokenTrack) break;
    }
    if (isBrokenTrack) return false;

    // All cuts passed, this is likely a cathode-crossing stopping muon.
    _isACathodeCrosser = true;
    return true;

  }

  // For MC events, check if the track is associated to a cosmic track
  bool StoppingMuonSelectionAlg::IsTrackMatchedToTrueCosmicTrack(art::Event const &evt,
                                                                 recob::PFParticle const &thisParticle) {
    const simb::MCParticle *particleP = 0x0;
    if (!evt.isRealData()) {
      particleP = truthUtil.GetMCParticleFromPFParticle(thisParticle,evt,fPFParticleTag);
      if (particleP == 0x0) return false;
      if ((pi_serv->TrackIdToMCTruth_P(particleP->TrackId())->Origin()) != simb::kCosmicRay)
        return true;
      else
        return false;
    }
    else
      return false;
  }

  // Order reco start and end point based on Y position
  void StoppingMuonSelectionAlg::OrderRecoStartEnd(TVector3 &start, TVector3 &end) {
    TVector3 prov;
    if (end.Y() > start.Y()) {
      prov = start;
      start = end;
      end = prov;
    }
    return;
  }

  // Get track from PFParticle
  const recob::Track StoppingMuonSelectionAlg::GetTrackFromPFParticle(art::Event const &evt,
                                                                      recob::PFParticle const &thisParticle) {
    const recob::Track *trackP = pfpUtil.GetPFParticleTrack(thisParticle,evt,fPFParticleTag,fTrackerTag);
    return *(trackP);
  }

  // Determine min and max hit peak time for this track
  void StoppingMuonSelectionAlg::SetMinAndMaxHitPeakTime(art::Event const &evt,
                                                         recob::PFParticle const &thisParticle,
                                                         double &minHitPeakTime,
                                                         double &maxHitPeakTime) {
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
