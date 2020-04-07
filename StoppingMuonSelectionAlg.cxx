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
    if (trackP == nullptr) {
      _isPFParticleATrack = false;
      return false;
    }
    _isPFParticleATrack = true;
    return true;
  }

  // Read parameters from FHICL file
  void StoppingMuonSelectionAlg::reconfigure(fhicl::ParameterSet const &p) {
    fTrackerTag = p.get<std::string>("TrackerTag", "pandoraTrack");
    fPFParticleTag = p.get<std::string>("PFParticleTag", "pandora");
    // Set cuts for CC.
    length_cutoff_CC               = p.get<double>("length_cutoff_CC",100);
    offsetFiducialBounds_CC        = p.get<double>("offsetFiducialBounds_CC",50);
    thicknessStartVolume_CC        = p.get<double>("thicknessStartVolume_CC",40);
    cutMinHitPeakTime_CC           = p.get<double>("cutMinHitPeakTime_CC",500);
    cutMaxHitPeakTime_CC           = p.get<double>("cutMaxHitPeakTime_CC",4800);
    radiusBrokenTracksSearch_CC    = p.get<double>("radiusBrokenTracksSearch_CC",50);
    cutCosAngleBrokenTracks_CC     = p.get<double>("cutCosAngleBrokenTracks_CC",0.995);
    cutCosAngleAlignment_CC        = p.get<double>("cutCosAngleAlignment_CC",0.995);
    cutContourAPA_CC               = p.get<double>("cutContourAPA_CC",10);
    // Set cuts for AC.
    length_cutoff_AC               = p.get<double>("length_cutoff_AC",100);
    offsetYStartPoint_AC           = p.get<double>("offsetYStartPoint_AC", 30);
    offsetZStartPoint_AC           = p.get<double>("offsetZStartPoint_AC", 50);
    cutMinHitPeakTime_AC           = p.get<double>("cutMinHitPeakTime_AC", 700);
    cutMaxHitPeakTime_AC           = p.get<double>("cutMaxHitPeakTime_AC", 4800);
    radiusBrokenTracksSearch_AC    = p.get<double>("radiusBrokenTracksSearch_AC", 10);
    cutCosAngleBrokenTracks_AC     = p.get<double>("cutCosAngleBrokenTracks_AC", 0.995);
    cutCosAngleAlignment_AC        = p.get<double>("cutCosAngleAlignment_AC", 0.995);
    cutContourAPA_AC               = p.get<double>("cutContourAPA_AC", 10);
    offsetFiducialBounds_AC        = p.get<double>("offsetFiducialBounds_AC",50);
    // Prepare geometry helper
    geoHelper.InitActiveVolumeBounds();
  }

  // Determine if the PFParticle is a selected cathode crosser
  bool StoppingMuonSelectionAlg::IsStoppingAnodeCrosser(art::Event const &evt,
                                                        recob::PFParticle const &thisParticle) {
    bool DEBUG = true;
    Reset();
    _evNumber = evt.id().event();

    // Set the correct boundaries for cathode crossers.
    geoHelper.SetFiducialBoundOffset(offsetFiducialBounds_AC);
    geoHelper.InitFiducialVolumeBounds();

    // Get recob::Track from PFParticle
    const recob::Track &track = GetTrackFromPFParticle(evt,thisParticle);
    _recoEndPoint = track.End<TVector3>();
    _recoStartPoint.SetXYZ(track.LocationAtPoint(track.FirstValidPoint()).X(), track.LocationAtPoint(track.FirstValidPoint()).Y(), track.LocationAtPoint(track.FirstValidPoint()).Z());
    OrderRecoStartEnd(_recoStartPoint, _recoEndPoint);
    _trackLength = track.Length();
    _trackID = track.ID();

    std::cout << "Track ID: " << _trackID << std::endl;

    // using the ordered start and end points calculate the angles _theta_xz and _theta_yz
    _theta_xz = TMath::RadToDeg() * TMath::ATan2(_recoStartPoint.X()-_recoEndPoint.X(), _recoStartPoint.Z()-_recoEndPoint.Z());
    _theta_yz = TMath::RadToDeg() * TMath::ATan2(_recoStartPoint.Y()-_recoEndPoint.Y(), _recoStartPoint.Z()-_recoEndPoint.Z());

    SetMinAndMaxHitPeakTime(evt,thisParticle,_minHitPeakTime,_maxHitPeakTime);

    if (_trackLength < length_cutoff_AC) {
      if (DEBUG) std::cout << "Track too short." << std::endl;
      return false;
    }

    if (TMath::Abs(_theta_yz-90)<10 || TMath::Abs(_theta_yz+90)<10 || TMath::Abs(_theta_xz-90)<10 || TMath::Abs(_theta_xz+90)<10) {
      if(DEBUG) std::cout << "Track's angle not accepted." << std::endl;
      return false;
    }

    bool goodTrack = geoHelper.IsPointYZProjectionInArea(_recoStartPoint,offsetYStartPoint_AC,offsetZStartPoint_AC);
    if (!goodTrack) {
      if(DEBUG) std::cout << "Track's start point not accepted." << std::endl;
      return false;
    }

    if (_minHitPeakTime <= cutMinHitPeakTime_AC) {
      if(DEBUG) std::cout << "Track's minHitPeakTime not accepted." << std::endl;
      return false;
    }
    if (_maxHitPeakTime >= cutMaxHitPeakTime_AC) {
      if(DEBUG) std::cout << "Track's maxHitPeakTime not accepted." << std::endl;
      return false;
    }

    if ((TMath::Abs(_recoStartPoint.Z()-geoHelper.GetAPABoundaries()[0])<=cutContourAPA_AC) ||
       (TMath::Abs(_recoStartPoint.Z()-geoHelper.GetAPABoundaries()[1])<=cutContourAPA_AC) ||
       (TMath::Abs(_recoEndPoint.Z()-geoHelper.GetAPABoundaries()[0])<=cutContourAPA_AC)   ||
       (TMath::Abs(_recoEndPoint.Z()-geoHelper.GetAPABoundaries()[1])<=cutContourAPA_AC)
     ) {
      if(DEBUG) std::cout << "Track's extreme points too close to APA." << std::endl;
      return false;
    }

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
      else {
        dirHigherTrack = dirSecondTrack;
        dirLowerTrack = dirFirstTrack;
        startPointLowerTrack = _recoStartPoint;
        endPointLowerTrack = _recoEndPoint;
        endPointHigherTrack = recoEndPointSecond;
      }

      TVector3 middlePointLowerTrack = (startPointLowerTrack + endPointLowerTrack) * 0.5;
      TVector3 dirHigherTrack_YZ(0., dirHigherTrack.Y(), dirHigherTrack.Z());
      TVector3 dirJoiningSegment(0., middlePointLowerTrack.Y()-endPointHigherTrack.Y(), middlePointLowerTrack.Z()-endPointHigherTrack.Z());
      double cosBeta = TMath::Cos(dirHigherTrack_YZ.Angle(dirJoiningSegment));
      double absCosAlpha = TMath::Abs(TMath::Cos(dirFirstTrack.Angle(dirSecondTrack)));

      if ( (absCosAlpha > cutCosAngleBrokenTracks_AC) && (cosBeta >= cutCosAngleAlignment_AC) )
        isBrokenTrack = true;

      double distHigherLower = TMath::Sqrt(TMath::Power(endPointHigherTrack.Y()-startPointLowerTrack.Y(),2) + TMath::Power(endPointHigherTrack.Z()-startPointLowerTrack.Z(),2));
      if (distHigherLower < radiusBrokenTracksSearch_AC) {
        if ( TMath::Abs(TMath::Cos(dirFirstTrack.Angle(dirSecondTrack))) > 0.96)
        isBrokenTrack = true;
      }
      // Stop searching if found one
      if (isBrokenTrack) break;
    }
    if (isBrokenTrack) {
      if(DEBUG) std::cout << "Track is broken." << std::endl;
      return false;
    }

    // Get the T0
    std::vector<anab::T0> pfparticleT0s = pfpUtil.GetPFParticleT0(thisParticle,evt,fPFParticleTag);
    if (pfparticleT0s.size() == 0) {
      _trackT0 = INV_DBL;
    }
    else
      _trackT0 = pfparticleT0s[0].Time();

    if (_trackT0 == INV_DBL) {
      if (CorrectPosAndGetT0(_recoStartPoint,_recoEndPoint) == INV_DBL) return false;
      CorrectPosEnd(_recoEndPoint, _minHitPeakTime, _maxHitPeakTime);
      trackInfo.isAnodeCrosserPandora = false;
      trackInfo.isAnodeCrosserMine = true;
    }
    else {
      std::cout << "Track tagged by Pandora." << std::endl;
      CorrectPosEnd(_recoEndPoint, _minHitPeakTime, _maxHitPeakTime);
      trackInfo.isAnodeCrosserPandora = true;
      trackInfo.isAnodeCrosserMine = false;
    }

    // Check if there are hits on the cryostat side if the track is selected by Pandora.
    if (trackInfo.isAnodeCrosserPandora) {
      // Get hit vector.
      auto const &trackHits = pfpUtil.GetPFParticleHits(thisParticle,evt,fPFParticleTag);
      if (!hitHelper.AreThereHitsOnCryoSide(trackHits)) {
        if(DEBUG) std::cout << "There are no hits in the cryostat sides." << std::endl;
        return false;
      }
    }

    bool goodEndPoint = geoHelper.IsPointInVolume(geoHelper.GetFiducialVolumeBounds(),_recoEndPoint);
    if (!goodEndPoint) {
      if(DEBUG) std::cout << "Track's end point not accepted." << std::endl;
      return false;
    }

    // All cuts passed, this is likely an anode-crossing stopping muon.
    _isAnAnodeCrosser = true;
    if (DEBUG) std::cout << "Track passed all selection cuts." << std::endl;
    return true;
  }

  // Correct position of the end point using the minimum hit peak time.
  void StoppingMuonSelectionAlg::CorrectPosEnd(TVector3 &_recoEndPoint, const double &minHitPeakTime, const double &maxHitPeakTime) {

    // To do: make this independent from protoDUNE run 1 setup.
    //double t0 = (minHitPeakTime - 500) / 2. * 1000.; // in ns.
    double drift_time_end_point = (maxHitPeakTime-minHitPeakTime) / 2. * 1000; // in ns.
    double drift_velocity = detprop->DriftVelocity()*1e-3; // in cm/ns.

    _recoEndPoint.SetX(drift_time_end_point * drift_velocity);

  }

  // Work out t0 for anode crossers.
  double StoppingMuonSelectionAlg::CorrectPosAndGetT0(TVector3 &_recoStartPoint, TVector3 &_recoEndPoint) {

    double drift_velocity = detprop->DriftVelocity()*1e-3; // in cm/ns.
    //std::cout << "StoppingMuonSelectionAlg::CorrectPosAndGetT0: " << "drift velocity = " << drift_velocity << std::endl;
    double drift_distance = geoHelper.GetAbsolutePlaneCoordinate(0); // First induction plane coordinate.
    //std::cout << "StoppingMuonSelectionAlg::CorrectPosAndGetT0: " << "drift distance = " << drift_distance << std::endl;

    if (_recoStartPoint.X() <= _recoEndPoint.X()) {
      if (_recoStartPoint.X() <= 0) {
        _trackT0 = (drift_distance - TMath::Abs(_recoStartPoint.X())) / drift_velocity;
        _recoStartPoint.SetX(_recoStartPoint.X() - (drift_velocity * _trackT0));
        _recoEndPoint.SetX(_recoEndPoint.X() - (drift_velocity * _trackT0));
      }
    }
    else {
      if (_recoStartPoint.X() >= 0) {
        _trackT0 = (drift_distance - TMath::Abs(_recoStartPoint.X())) / drift_velocity;
        _recoStartPoint.SetX(_recoStartPoint.X() + (drift_velocity * _trackT0));
        _recoEndPoint.SetX(_recoEndPoint.X() + (drift_velocity * _trackT0));
      }
    }

    return _trackT0;
  }

  // Determine if the PFParticle is a selected cathode crosser
  bool StoppingMuonSelectionAlg::IsStoppingCathodeCrosser(art::Event const &evt,
                                                          recob::PFParticle const &thisParticle) {
    Reset();
    _evNumber = evt.id().event();

    // Set the correct boundaries for cathode crossers.
    geoHelper.SetFiducialBoundOffset(offsetFiducialBounds_CC);
    geoHelper.SetThicknessStartVolume(thicknessStartVolume_CC);
    geoHelper.InitFiducialVolumeBounds();

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
    //std::cout << "theta xz :"  << _theta_xz << std::endl;
    //std::cout << "theta yz :"  << _theta_yz << std::endl;
    SetMinAndMaxHitPeakTime(evt,thisParticle,_minHitPeakTime,_maxHitPeakTime);
    //std::cout << "_minHitPeakTime :"  << _minHitPeakTime << std::endl;
    //std::cout << "_maxHitPeakTime :"  << _maxHitPeakTime << std::endl;
    // Apply cuts with selection with progressive cuts
    //std::cout << "Length: " << _trackLength << std::endl;

    // Get the T0
    std::vector<anab::T0> pfparticleT0s = pfpUtil.GetPFParticleT0(thisParticle,evt,fPFParticleTag);
    if (pfparticleT0s.size() == 0) {
      _trackT0 = INV_DBL;
      return false;
    }
    else
      _trackT0 = pfparticleT0s[0].Time();

    if (_trackLength < length_cutoff_CC) return false;
    if (TMath::Abs(_theta_yz-90)<10 || TMath::Abs(_theta_yz+90)<10 || TMath::Abs(_theta_xz-90)<10 || TMath::Abs(_theta_xz+90)<10) return false;
    bool goodTrack = ((_recoStartPoint.X()*_recoEndPoint.X()<0)
                     && geoHelper.IsPointInSlice(_recoStartPoint)
                     && geoHelper.IsPointInVolume(geoHelper.GetFiducialVolumeBounds(),_recoEndPoint));
    //std::cout << "start*end: " << _recoStartPoint.X()*_recoEndPoint.X() << std::endl;
    //std::cout << "start point in slice: " << geoHelper.IsPointInSlice(_recoStartPoint) << std::endl;
    //std::cout << "End point in volume: " << geoHelper.IsPointInVolume(geoHelper.GetFiducialVolumeBounds(),_recoEndPoint) << std::endl;
    //std::cout << "End point (XYZ): " << _recoEndPoint.X() << " " << _recoEndPoint.Y() << " " << _recoEndPoint.Z() << std::endl;
    if (!goodTrack) return false;
    if (_minHitPeakTime <= cutMinHitPeakTime_CC) return false;
    if (_maxHitPeakTime >= cutMaxHitPeakTime_CC) return false;
    if ((TMath::Abs(_recoEndPoint.Z()-geoHelper.GetAPABoundaries()[0])<=cutContourAPA_CC) || (TMath::Abs(_recoEndPoint.Z()-geoHelper.GetAPABoundaries()[1])<=cutContourAPA_CC)) return false;
    //std::cout << "End Point Z: " << _recoEndPoint.Z() << std::endl;
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
      TVector3 middlePointLowerTrack = (startPointLowerTrack + endPointLowerTrack) * 0.5;
      TVector3 dirHigherTrack_YZ(0., dirHigherTrack.Y(), dirHigherTrack.Z());
      TVector3 dirJoiningSegment(0., middlePointLowerTrack.Y()-endPointHigherTrack.Y(), middlePointLowerTrack.Z()-endPointHigherTrack.Z());
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
    trackInfo.isCathodeCrosser = true;
    return true;

  }

  // For MC events, check if the track is associated to a cosmic track
  bool StoppingMuonSelectionAlg::IsTrackMatchedToTrueCosmicTrack(art::Event const &evt,
                                                                 recob::PFParticle const &thisParticle) {
    const simb::MCParticle *particleP = 0x0;
    if (!evt.isRealData()) {
      particleP = truthUtil.GetMCParticleFromPFParticle(thisParticle,evt,fPFParticleTag);
      if (particleP == 0x0) return false;
      if ((pi_serv->TrackIdToMCTruth_P(particleP->TrackId())->Origin()) == simb::kCosmicRay) {
        //std::cout << "Particle ID: " << particleP->TrackId() << " Pdg: " << particleP->PdgCode() << std::endl;
        //int pdg_mother = pi_serv->TrackIdToMotherParticle_P(particleP->TrackId())->PdgCode();
        //std::cout << "StoppingMuonSelectionAlg::IsTrackMatchedToTrueCosmicTrack: " << "pdg of mother particle: " << pdg_mother << std::endl;
        //std::cout << "process: " << particleP->Process() << std::endl;
        return true;
      }
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

  // Set MCParticle properties
  void StoppingMuonSelectionAlg::SetMCParticleProperties(art::Event const &evt, recob::PFParticle const &thisParticle) {
    const simb::MCParticle *particleP = truthUtil.GetMCParticleFromPFParticle(thisParticle,evt,fPFParticleTag);
    _pdg = particleP->PdgCode();
    double *av = geoHelper.GetActiveVolumeBounds();
    int firstPoint = truthUtil.GetFirstTrajectoryPointInTPCActiveVolume(*particleP,av[0],av[1],av[2],av[3],av[4],av[5]);
    _trueStartPoint.SetXYZ(particleP->Vx(firstPoint),particleP->Vy(firstPoint),particleP->Vz(firstPoint));
    _trueEndPoint = particleP->EndPosition().Vect();
    _trueStartT = particleP->T(firstPoint);
    _trueEndT = particleP->EndPosition().T();
    _trueTrackID = particleP->TrackId();
    _areMCParticlePropertiesSet = true;
  }

  // Check if the true track associated with that PFParticle is a stopping muon
  bool StoppingMuonSelectionAlg::IsTrueParticleACathodeCrossingStoppingMuon(art::Event const &evt, recob::PFParticle const &thisParticle) {
    if (!_areMCParticlePropertiesSet)
      SetMCParticleProperties(evt,thisParticle);
    return (TMath::Abs(_pdg)==13
            && (_trueStartPoint.X()*_trueEndPoint.X()<0)
            && geoHelper.IsPointInVolume(geoHelper.GetActiveVolumeBounds(),_trueEndPoint));
  }

  // Check if the true track associated with that PFParticle is a stopping muon
  bool StoppingMuonSelectionAlg::IsTrueParticleAnAnodeCrossingStoppingMuon(art::Event const &evt, recob::PFParticle const &thisParticle) {
    if (!_areMCParticlePropertiesSet)
      SetMCParticleProperties(evt,thisParticle);
    return (TMath::Abs(_pdg)==13 &&
           (TMath::Abs(_trueStartPoint.X()-geoHelper.GetActiveVolumeBounds()[1]) < 1. || TMath::Abs(_trueStartPoint.X()+geoHelper.GetActiveVolumeBounds()[1]<1.)) &&
           geoHelper.IsPointInVolume(geoHelper.GetActiveVolumeBounds(),_trueEndPoint));
  }

  // N-1 cuts for Cathode crossers
  bool StoppingMuonSelectionAlg::NMinus1Cathode(const std::string &excludeCut, art::Event const &evt, const recob::PFParticle &thisParticle) {

    Reset();
    _evNumber = evt.id().event();

    // Set the correct boundaries for cathode crossers.
    geoHelper.SetFiducialBoundOffset(offsetFiducialBounds_CC);
    geoHelper.SetThicknessStartVolume(thicknessStartVolume_CC);
    geoHelper.InitFiducialVolumeBounds();

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
    // Apply cuts with selection with progressive cuts
    if (_trackLength < length_cutoff_CC) return false;

    // Apply cuts for good track separately.
    if (!(_recoStartPoint.X()*_recoEndPoint.X()<0)) return false;

    if (excludeCut!="thicknessStartVolume") {
      bool goodStart = geoHelper.IsPointInSlice(_recoStartPoint);
      if (!goodStart) return false;
    }
    double *fidBounds = geoHelper.GetFiducialVolumeBounds();
    if (excludeCut != "distanceFiducialVolumeX" && excludeCut != "distanceFiducialVolumeY" && excludeCut != "distanceFiducialVolumeZ")   {
      bool goodEndPoint = geoHelper.IsPointInVolume(geoHelper.GetFiducialVolumeBounds(), _recoEndPoint);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeX") {
      bool goodEndPoint = (_recoEndPoint.Y() >= fidBounds[2] && _recoEndPoint.Y() <= fidBounds[3] && _recoEndPoint.Z() >= fidBounds[4] && _recoEndPoint.Z() <= fidBounds[5]);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeY") {
      bool goodEndPoint = (_recoEndPoint.X() >= fidBounds[0] && _recoEndPoint.X() <= fidBounds[1] && _recoEndPoint.Z() >= fidBounds[4] && _recoEndPoint.Z() <= fidBounds[5]);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeZ") {
      bool goodEndPoint = (_recoEndPoint.X() >= fidBounds[0] && _recoEndPoint.X() <= fidBounds[1] && _recoEndPoint.Y() >= fidBounds[2] && _recoEndPoint.Y() <= fidBounds[3]);
      if (!goodEndPoint) return false;
    }

    if (excludeCut!="cutMinHitPeakTime") {
      if (_minHitPeakTime <= cutMinHitPeakTime_CC) return false;
    }
    if (excludeCut!="cutMaxHitPeakTime") {
      if (_maxHitPeakTime >= cutMaxHitPeakTime_CC) return false;
    }

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
      TVector3 middlePointLowerTrack = (startPointLowerTrack + endPointLowerTrack) * 0.5;
      TVector3 dirHigherTrack_YZ(0., dirHigherTrack.Y(), dirHigherTrack.Z());
      TVector3 dirJoiningSegment(0., middlePointLowerTrack.Y()-endPointHigherTrack.Y(), middlePointLowerTrack.Z()-endPointHigherTrack.Z());
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
    trackInfo.isCathodeCrosser = true;
    return true;
  }

  // N-1 cuts for Anode crossers
  bool StoppingMuonSelectionAlg::NMinus1Anode(const std::string &excludeCut, art::Event const &evt, const recob::PFParticle &thisParticle) {

    bool DEBUG = false;
    Reset();
    _evNumber = evt.id().event();

    // Set the correct boundaries for cathode crossers.
    geoHelper.SetFiducialBoundOffset(offsetFiducialBounds_AC);
    geoHelper.InitFiducialVolumeBounds();

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

    SetMinAndMaxHitPeakTime(evt,thisParticle,_minHitPeakTime,_maxHitPeakTime);

    if (_trackLength < length_cutoff_AC) {
      if (DEBUG) std::cout << "Track too short." << std::endl;
      return false;
    }

    if (TMath::Abs(_theta_yz-90)<10 || TMath::Abs(_theta_yz+90)<10 || TMath::Abs(_theta_xz-90)<10 || TMath::Abs(_theta_xz+90)<10) {
      if(DEBUG) std::cout << "Track's angle not accepted." << std::endl;
      return false;
    }

    bool goodTrack = false;
    double *actBounds = geoHelper.GetActiveVolumeBounds();
    if (excludeCut != "offsetYStartPoint" && excludeCut != "offsetZStartPoint")
      goodTrack = geoHelper.IsPointYZProjectionInArea(_recoStartPoint,offsetYStartPoint_AC,offsetZStartPoint_AC);
    else if (excludeCut == "offsetYStartPoint")
      goodTrack = (_recoStartPoint.Z()>=(actBounds[4]+offsetZStartPoint_AC) && _recoStartPoint.Z()<=(actBounds[5]-offsetZStartPoint_AC));
    else if (excludeCut == "offsetZStartPoint")
      goodTrack = (_recoStartPoint.Y()>=(actBounds[2]+offsetYStartPoint_AC) && _recoStartPoint.Y()<=(actBounds[3]-offsetYStartPoint_AC));
    if (!goodTrack) {
      if(DEBUG) std::cout << "Track's start point not accepted." << std::endl;
      return false;
    }

    if (excludeCut != "cutMinHitPeakTime") {
      if (_minHitPeakTime <= cutMinHitPeakTime_AC) {
        if(DEBUG) std::cout << "Track's minHitPeakTime not accepted." << std::endl;
        return false;
      }
    }

    if (excludeCut != "cutMaxHitPeakTime") {
      if (_maxHitPeakTime >= cutMaxHitPeakTime_AC) {
        if(DEBUG) std::cout << "Track's maxHitPeakTime not accepted." << std::endl;
        return false;
      }
    }

    if ((TMath::Abs(_recoStartPoint.Z()-geoHelper.GetAPABoundaries()[0])<=cutContourAPA_AC) ||
       (TMath::Abs(_recoStartPoint.Z()-geoHelper.GetAPABoundaries()[1])<=cutContourAPA_AC) ||
       (TMath::Abs(_recoEndPoint.Z()-geoHelper.GetAPABoundaries()[0])<=cutContourAPA_AC)   ||
       (TMath::Abs(_recoEndPoint.Z()-geoHelper.GetAPABoundaries()[1])<=cutContourAPA_AC)
     ) {
      if(DEBUG) std::cout << "Track's extreme points too close to APA." << std::endl;
      return false;
    }

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
      else {
        dirHigherTrack = dirSecondTrack;
        dirLowerTrack = dirFirstTrack;
        startPointLowerTrack = _recoStartPoint;
        endPointLowerTrack = _recoEndPoint;
        endPointHigherTrack = recoEndPointSecond;
      }

      TVector3 middlePointLowerTrack = (startPointLowerTrack + endPointLowerTrack) * 0.5;
      TVector3 dirHigherTrack_YZ(0., dirHigherTrack.Y(), dirHigherTrack.Z());
      TVector3 dirJoiningSegment(0., middlePointLowerTrack.Y()-endPointHigherTrack.Y(), middlePointLowerTrack.Z()-endPointHigherTrack.Z());
      double cosBeta = TMath::Cos(dirHigherTrack_YZ.Angle(dirJoiningSegment));
      double absCosAlpha = TMath::Abs(TMath::Cos(dirFirstTrack.Angle(dirSecondTrack)));

      if ( (absCosAlpha > cutCosAngleBrokenTracks_AC) && (cosBeta >= cutCosAngleAlignment_AC) )
        isBrokenTrack = true;

      double distHigherLower = TMath::Sqrt(TMath::Power(endPointHigherTrack.Y()-startPointLowerTrack.Y(),2) + TMath::Power(endPointHigherTrack.Z()-startPointLowerTrack.Z(),2));
      if (distHigherLower < radiusBrokenTracksSearch_AC) {
        if ( TMath::Abs(TMath::Cos(dirFirstTrack.Angle(dirSecondTrack))) > 0.96)
        isBrokenTrack = true;
      }
      // Stop searching if found one
      if (isBrokenTrack) break;
    }
    if (isBrokenTrack) {
      if(DEBUG) std::cout << "Track is broken." << std::endl;
      return false;
    }

    // Get the T0
    std::vector<anab::T0> pfparticleT0s = pfpUtil.GetPFParticleT0(thisParticle,evt,fPFParticleTag);
    if (pfparticleT0s.size() == 0) {
      _trackT0 = INV_DBL;
    }
    else
      _trackT0 = pfparticleT0s[0].Time();

    if (_trackT0 == INV_DBL) {
      if (CorrectPosAndGetT0(_recoStartPoint,_recoEndPoint) == INV_DBL) return false;
      //CorrectPosEnd(_recoStartPoint, _recoEndPoint, _minHitPeakTime);
      trackInfo.isAnodeCrosserPandora = false;
      trackInfo.isAnodeCrosserMine = true;
    }
    else {
      if (DEBUG) std::cout << "Track tagged by Pandora." << std::endl;
      //CorrectPosEnd(_recoStartPoint, _recoEndPoint, _minHitPeakTime);
      trackInfo.isAnodeCrosserPandora = true;
      trackInfo.isAnodeCrosserMine = false;
    }

    // Check if there are hits on the cryostat side if the track is selected by Pandora.
    if (trackInfo.isAnodeCrosserPandora) {
      // Get hit vector.
      auto const &trackHits = pfpUtil.GetPFParticleHits(thisParticle,evt,fPFParticleTag);
      if (!hitHelper.AreThereHitsOnCryoSide(trackHits)) {
        if(DEBUG) std::cout << "There are no hits in the cryostat sides." << std::endl;
        return false;
      }
    }

    double *fidBounds = geoHelper.GetFiducialVolumeBounds();
    if (excludeCut != "distanceFiducialVolumeX" && excludeCut != "distanceFiducialVolumeY" && excludeCut != "distanceFiducialVolumeZ")   {
      bool goodEndPoint = geoHelper.IsPointInVolume(geoHelper.GetFiducialVolumeBounds(), _recoEndPoint);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeX") {
      bool goodEndPoint = (_recoEndPoint.Y() >= fidBounds[2] && _recoEndPoint.Y() <= fidBounds[3] && _recoEndPoint.Z() >= fidBounds[4] && _recoEndPoint.Z() <= fidBounds[5]);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeY") {
      bool goodEndPoint = (_recoEndPoint.X() >= fidBounds[0] && _recoEndPoint.X() <= fidBounds[1] && _recoEndPoint.Z() >= fidBounds[4] && _recoEndPoint.Z() <= fidBounds[5]);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeZ") {
      bool goodEndPoint = (_recoEndPoint.X() >= fidBounds[0] && _recoEndPoint.X() <= fidBounds[1] && _recoEndPoint.Y() >= fidBounds[2] && _recoEndPoint.Y() <= fidBounds[3]);
      if (!goodEndPoint) return false;
    }

    // All cuts passed, this is likely an anode-crossing stopping muon.
    _isAnAnodeCrosser = true;
    if (DEBUG) std::cout << "Track passed all selection cuts." << std::endl;
    return true;

  }

  // N-1 cuts for Cathode crossers, simple version
  bool StoppingMuonSelectionAlg::NMinus1CathodeSimple(const std::string &excludeCut, art::Event const &evt, const recob::PFParticle &thisParticle) {

    Reset();
    _evNumber = evt.id().event();

    // Set the correct boundaries for cathode crossers.
    geoHelper.SetFiducialBoundOffset(offsetFiducialBounds_CC);
    geoHelper.SetThicknessStartVolume(thicknessStartVolume_CC);
    geoHelper.InitFiducialVolumeBounds();

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
    // Apply cuts with selection with progressive cuts
    if (_trackLength < length_cutoff_CC) return false;

    // Apply cuts for good track separately.
    if (!(_recoStartPoint.X()*_recoEndPoint.X()<0)) return false;

    if (excludeCut!="thicknessStartVolume") {
      bool goodStart = geoHelper.IsPointInSlice(_recoStartPoint);
      if (!goodStart) return false;
    }
    double *fidBounds = geoHelper.GetFiducialVolumeBounds();
    if (excludeCut != "distanceFiducialVolumeX" && excludeCut != "distanceFiducialVolumeY" && excludeCut != "distanceFiducialVolumeZ")   {
      bool goodEndPoint = geoHelper.IsPointInVolume(geoHelper.GetFiducialVolumeBounds(), _recoEndPoint);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeX") {
      bool goodEndPoint = (_recoEndPoint.Y() >= fidBounds[2] && _recoEndPoint.Y() <= fidBounds[3] && _recoEndPoint.Z() >= fidBounds[4] && _recoEndPoint.Z() <= fidBounds[5]);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeY") {
      bool goodEndPoint = (_recoEndPoint.X() >= fidBounds[0] && _recoEndPoint.X() <= fidBounds[1] && _recoEndPoint.Z() >= fidBounds[4] && _recoEndPoint.Z() <= fidBounds[5]);
      if (!goodEndPoint) return false;
    }
    else if (excludeCut == "distanceFiducialVolumeZ") {
      bool goodEndPoint = (_recoEndPoint.X() >= fidBounds[0] && _recoEndPoint.X() <= fidBounds[1] && _recoEndPoint.Y() >= fidBounds[2] && _recoEndPoint.Y() <= fidBounds[3]);
      if (!goodEndPoint) return false;
    }

    if (excludeCut!="cutMinHitPeakTime") {
      if (_minHitPeakTime <= cutMinHitPeakTime_CC) return false;
    }
    if (excludeCut!="cutMaxHitPeakTime") {
      if (_maxHitPeakTime >= cutMaxHitPeakTime_CC) return false;
    }

    // All cuts passed, this is likely a cathode-crossing stopping muon.
    _isACathodeCrosser = true;
    trackInfo.isCathodeCrosser = true;
    return true;
  }

  // Get the property for this track. Only if its selected.
  const trackProperties StoppingMuonSelectionAlg::GetTrackProperties() {
    trackInfo.evNumber = _evNumber;
    trackInfo.trackT0 = _trackT0;
    trackInfo.recoStartPoint = _recoStartPoint;
    trackInfo.recoEndPoint = _recoEndPoint;
    trackInfo.theta_xz = _theta_xz;
    trackInfo.theta_yz = _theta_yz;
    trackInfo.minHitPeakTime = _minHitPeakTime;
    trackInfo.maxHitPeakTime = _maxHitPeakTime;
    trackInfo.trackLength = _trackLength;
    trackInfo.trackID = _trackID;
    trackInfo.pdg = _pdg;
    trackInfo.trueStartPoint = _trueStartPoint;
    trackInfo.trueEndPoint = _trueEndPoint;
    trackInfo.trueStartT = _trueStartT;
    trackInfo.trueEndT = _trueEndT;
    trackInfo.trueTrackID = _trueTrackID;
    return trackInfo;
  }

  void StoppingMuonSelectionAlg::Reset() {
    _isACathodeCrosser = false;
    _isAnAnodeCrosser = false;
    _isPFParticleATrack = false;
    _areMCParticlePropertiesSet = false;
    _evNumber = INV_INT;
    _trackT0 = INV_DBL;
    _recoStartPoint.SetXYZ(INV_DBL,INV_DBL,INV_DBL);
    _recoEndPoint.SetXYZ(INV_DBL,INV_DBL,INV_DBL);
    _theta_xz = INV_DBL;
    _theta_yz = INV_DBL;
    _minHitPeakTime = INV_DBL;
    _maxHitPeakTime = INV_DBL;
    _trackLength = INV_DBL;
    _trackID = INV_DBL;
    _pdg = INV_INT;
    _trueStartPoint.SetXYZ(INV_DBL,INV_DBL,INV_DBL);
    _trueEndPoint.SetXYZ(INV_DBL,INV_DBL,INV_DBL);
    _trueStartT = INV_DBL;
    _trueEndT = INV_DBL;
    _trueTrackID = INV_DBL;
    trackInfo.Reset();
  }

} // end of namespace stoppingcosmicmuonselection

#endif
