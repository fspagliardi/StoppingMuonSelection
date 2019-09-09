/***
  Class containing useful functions for geometry.

*/
#ifndef CALORIMETRY_HELPER_CXX
#define CALORIMETRY_HELPER_CXX

#include "CalorimetryHelper.h"

namespace stoppingcosmicmuonselection {

  CalorimetryHelper::CalorimetryHelper() {

  }

  CalorimetryHelper::CalorimetryHelper(const recob::PFParticle &thisParticle, art::Event const &evt) {
    Set(thisParticle, evt);
  }

  CalorimetryHelper::~CalorimetryHelper() {

  }

  // Get the calorimetry from the PFParticle
  void CalorimetryHelper::Set(const recob::PFParticle &thisParticle, art::Event const &evt) {
    Reset();
    // Set variable to see if it's data or MC (different histo scales)
    if (evt.isRealData()) _isData = true;
    const recob::Track &track = *(pfpUtil.GetPFParticleTrack(thisParticle,evt,fPFParticleTag,fTrackerTag));
    _calos = trackUtil.GetRecoTrackCalorimetry(track,evt,fTrackerTag,fCalorimetryTag);
    _isCalorimetrySet = true;
    if (!IsValid()) return;
    for (size_t itcal = 0; itcal < _calos.size(); itcal++) {
      if (!(_calos[itcal].PlaneID().isValid)) {std::cout << "plane not valid"<< std::endl;continue;}
      int planeNumb = _calos[itcal].PlaneID().Plane;
      if (planeNumb<0 || planeNumb>2) {std::cout << "plane number not valid"<< std::endl;continue;}
      size_t const Nhits = _calos[itcal].dEdx().size();
      _trackHitNumb[planeNumb] = int(Nhits);
      //geo::PlaneID Plane = _calos[itcal].PlaneID();
      for (size_t itHit = 0; itHit < Nhits; itHit++)  {
        //std::cout << "TpIndex: " << (_calos[itcal].TpIndices())[itHit] << std::endl;
        auto const & TrackPos = (_calos[itcal].XYZ())[itHit];
        _dqdx[planeNumb][itHit]=(_calos[itcal].dQdx())[itHit];
        _dedx[planeNumb][itHit]=(_calos[itcal].dEdx())[itHit];
        _resrange[planeNumb][itHit]=(_calos[itcal].ResidualRange())[itHit];
        _hitx[planeNumb][itHit]=TrackPos.X();
    	  _hity[planeNumb][itHit]=TrackPos.Y();
    	  _hitz[planeNumb][itHit]=TrackPos.Z();
        _track_pitch[planeNumb][itHit]=(_calos[itcal].TrkPitchVec())[itHit];
        // Get Hit Peak Time
        //const size_t & hitIndex = (_calos[itcal].TpIndices())[itHit];
        //const auto & thisHit = allHits[hitIndex];
        //hitPeakTime[planeNumb][itHit] = thisHit.PeakTime();
        _hitPeakTime[planeNumb][itHit] = INV_DBL; // For now in MCC11
        // Apply Lifetime corrections
        const geo::Point_t HitPoint(TrackPos.X(), TrackPos.Y(), TrackPos.Z());
        geo::TPCID const & tpcid = geom->FindTPCAtPosition(HitPoint);
        if (!tpcid.isValid) {std::cout << "tpc not valid"<< std::endl;continue;}
        int CryoID = geom->FindCryostatAtPosition(HitPoint);
        double Ticks = detprop->ConvertXToTicks(TrackPos.X(), planeNumb, tpcid.TPC, CryoID);
        _drift_time[planeNumb][itHit] = (Ticks - detprop->TriggerOffset()) * detprop->SamplingRate()*1e-3;
        _corr_factors[planeNumb][itHit] = LifeTimeCorr(Ticks, 0);
      }
    }
    OrderResRange();
  }

  // Check if the calorimetry is valid
  bool CalorimetryHelper::IsValid() {
    if (_calos.size() == 0) {
      _isValid = false;
      return false;
    }
    else {
      _isValid = true;
      return true;
    }
  }

  // Order the residual range with respect to the track direction
  void CalorimetryHelper::OrderResRange() {
    for (int planeNumber = 0; planeNumber < 3; planeNumber++) {
      std::vector<double> res_vect;
      for (int i = 0; i < TMath::Min(3000, _trackHitNumb[planeNumber]);i++) {
        res_vect.push_back(_resrange[planeNumber][i]);
      }
      int size = res_vect.size();
      double max = *max_element(res_vect.begin(),res_vect.end());
      for (int i = 0; i < TMath::Min(3000, _trackHitNumb[planeNumber]);i++) {
        if (_hity[planeNumber][size-1] < _hity[planeNumber][0])  {
          if (_resrange[planeNumber][size-1] < _resrange[planeNumber][0])  {
            _resrange_ord[planeNumber][i] = _resrange[planeNumber][i];
          }
          else {
            _resrange_ord[planeNumber][i] = max-_resrange[planeNumber][i];
          }
        }
        else {
          if (_resrange[planeNumber][size-1] < _resrange[planeNumber][0])  {
            _resrange_ord[planeNumber][i] = max-_resrange[planeNumber][i];
          }
          else {
            _resrange_ord[planeNumber][i] = _resrange[planeNumber][i];
          }
        }
      }
    } // end loop on planes
  }

  // Get hit numb
  int CalorimetryHelper::GetHitNumb(const int &planeNumb) {
    return _trackHitNumb[planeNumb];
  }
  // Get dqdx
  double *CalorimetryHelper::GetdQdx(const int &planeNumb) {
    return _dqdx[planeNumb];
  }
  // Get dEdx
  double *CalorimetryHelper::GetdEdx(const int &planeNumb) {
    return _dedx[planeNumb];
  }
  // Get Residual range
  double *CalorimetryHelper::GetResRange(const int &planeNumb) {
    return _resrange[planeNumb];
  }
  // Get Ordered residual range
  double *CalorimetryHelper::GetResRangeOrdered(const int &planeNumb) {
    return _resrange_ord[planeNumb];
  }
  // Get HitX
  double *CalorimetryHelper::GetHitX(const int &planeNumb) {
    return _hitx[planeNumb];
  }
  // Get HitY
  double *CalorimetryHelper::GetHitY(const int &planeNumb) {
    return _hity[planeNumb];
  }
  // Get HitZ
  double *CalorimetryHelper::GetHitZ(const int &planeNumb) {
    return _hitz[planeNumb];
  }
  // Get HitPeakTime
  double *CalorimetryHelper::GetHitPeakTime(const int &planeNumb) {
    return _hitPeakTime[planeNumb];
  }
  // Get lifetime correction factors
  double *CalorimetryHelper::GetCorrFactor(const int &planeNumb) {
    return _corr_factors[planeNumb];
  }
  // Get drift times
  double *CalorimetryHelper::GetDriftTime(const int &planeNumb) {
    return _drift_time[planeNumb];
  }
  // Get track pitches
  double *CalorimetryHelper::GetTrackPitch(const int &planeNumb) {
    return _track_pitch[planeNumb];
  }

  // Get the lifetime correction
  double CalorimetryHelper::LifeTimeCorr(double &ticks, const double &T0) {
    double timetick = detprop->SamplingRate()*1e-3; // Sample in microsec
    double presampling = detprop->TriggerOffset();
    //std::cout << "Presampling: " << presampling << std::endl;
    ticks = ticks - presampling;
    double time;
    time = ticks * timetick - T0;
    double tau = detprop->ElectronLifetime();
    double correction = TMath::Exp(time/tau);
    return correction;
  }

  // FIll 2D histo of dQdx vs residual range for hits in a given plane
  void CalorimetryHelper::FillHisto_dQdxVsRR(TH2D *h_dQdxVsRR, const int &planeNumb) {
    int hitNumb = GetHitNumb(planeNumb);
    double *dQdx = GetdQdx(planeNumb);
    double *resRangeOrd = GetResRangeOrdered(planeNumb);
    for (int it = 0; it < hitNumb; it++) {
      h_dQdxVsRR->Fill(resRangeOrd[it],dQdx[it]);
    }
  }

  // FIll 2D histo of dQdx vs residual range for hits in a given plane, in a track pitch interval
  void CalorimetryHelper::FillHisto_dQdxVsRR(TH2D *h_dQdxVsRR, const int &planeNumb, const double &tp_min, const double &tp_max) {
    int hitNumb = GetHitNumb(planeNumb);
    double *dQdx = GetdQdx(planeNumb);
    double *resRangeOrd = GetResRangeOrdered(planeNumb);
    double *trackPitch = GetTrackPitch(planeNumb);
    for (int it = 0; it < hitNumb; it++) {
      if (trackPitch[it]<tp_min || trackPitch[it]>tp_max) continue;
      h_dQdxVsRR->Fill(resRangeOrd[it],dQdx[it]);
    }
  }

  // FIll 2D histo of dQdx vs residual range for hits in a given plane. Correct by MC lifetime
  void CalorimetryHelper::FillHisto_dQdxVsRR_LTCorr(TH2D *h_dQdxVsRR, const int &planeNumb) {
    int hitNumb = GetHitNumb(planeNumb);
    double *dQdx = GetdQdx(planeNumb);
    double *resRangeOrd = GetResRangeOrdered(planeNumb);
    double *corrFactor = GetCorrFactor(planeNumb);
    for (int it = 0; it < hitNumb; it++)
      h_dQdxVsRR->Fill(resRangeOrd[it],dQdx[it]*corrFactor[it]);
  }

  // Same as above but with track pitch cut
  void CalorimetryHelper::FillHisto_dQdxVsRR_LTCorr(TH2D *h_dQdxVsRR, const int &planeNumb, const double &tp_min, const double &tp_max) {
    int hitNumb = GetHitNumb(planeNumb);
    double *dQdx = GetdQdx(planeNumb);
    double *resRangeOrd = GetResRangeOrdered(planeNumb);
    double *trackPitch = GetTrackPitch(planeNumb);
    double *corrFactor = GetCorrFactor(planeNumb);
    for (int it = 0; it < hitNumb; it++) {
      if (trackPitch[it]<tp_min || trackPitch[it]>tp_max) continue;
      h_dQdxVsRR->Fill(resRangeOrd[it],dQdx[it]*corrFactor[it]);
    }
  }

  // Set the parameters from the FHICL file
  void CalorimetryHelper::reconfigure(fhicl::ParameterSet const &p) {
    fTrackerTag = p.get<std::string>("TrackerTag");
    fCalorimetryTag = p.get<std::string>("CalorimetryTag");
    fPFParticleTag = p.get<std::string>("PFParticleTag");
  }

  // Reset
  void CalorimetryHelper::Reset() {
    _isValid = false;
    _isCalorimetrySet = false;
    _isData = false;
    _calos.clear();
    for (int j=0; j<3; j++) {
      _trackHitNumb[j]=INV_DBL;
      for(int k=0; k<3000; k++){
  	    _dqdx[j][k]=INV_DBL;
  	    _dedx[j][k]=INV_DBL;
  	    _resrange[j][k]=INV_DBL;
        _resrange_ord[j][k]=INV_DBL;
  	    _hitx[j][k]=INV_DBL;
  	    _hity[j][k]=INV_DBL;
  	    _hitz[j][k]=INV_DBL;
        _hitPeakTime[j][k]=INV_DBL;
        _corr_factors[j][k] = INV_DBL;
        _drift_time[j][k] = INV_DBL;
        _track_pitch[j][k] = INV_DBL;
  	  }
    }
  }

} // end of namespace stoppingcosmicmuonselection

#endif
