/***
  Class containing useful functions for calorimetry.

*/
#ifndef CALORIMETRY_HELPER_CXX
#define CALORIMETRY_HELPER_CXX

#include "CalorimetryHelper.h"

namespace stoppingcosmicmuonselection {

  CalorimetryHelper::CalorimetryHelper() {

  }

  CalorimetryHelper::CalorimetryHelper(const recob::PFParticle &thisParticle, art::Event const &evt, const int &plane) {
    Set(thisParticle, evt, plane);
  }

  CalorimetryHelper::~CalorimetryHelper() {

  }

  // Get the calorimetry from the PFParticle
  void CalorimetryHelper::Set(const recob::PFParticle &thisParticle, art::Event const &evt, const int &plane) {
    Reset();
    _plane = plane;
    bool correct_dQdx = false;
    // Set variable to see if it's data or MC (different histo scales)
    if (evt.isRealData()) _isData = true;
    const recob::Track &track = *(pfpUtil.GetPFParticleTrack(thisParticle,evt,fPFParticleTag,fTrackerTag));
    _calos = trackUtil.GetRecoTrackCalorimetry(track,evt,fTrackerTag,fCalorimetryTag);
    _isCalorimetrySet = true;

    if (!IsValid()) {
      std::cout << "CalorimetryHelper.cxx: Calorimetry invalid for plane: " << plane << std::endl;
      return;
    }

    for (size_t itcal = 0; itcal < _calos.size(); itcal++) {
      if (!(_calos[itcal].PlaneID().isValid)) {
        std::cout << "CalorimetryHelper.cxx: " << "plane at entry " << itcal << " not valid"<< std::endl;
        continue;
      }

      int planeNumb = _calos[itcal].PlaneID().Plane;
      if (plane != planeNumb) continue;
      size_t const Nhits = _calos[itcal].dEdx().size();
      _trackHitNumb = int(Nhits);
      //geo::PlaneID Plane = _calos[itcal].PlaneID();
      for (size_t itHit = 0; itHit < Nhits; itHit++)  {
        //std::cout << "TpIndex: " << (_calos[itcal].TpIndices())[itHit] << std::endl;
        auto const & TrackPos = (_calos[itcal].XYZ())[itHit];
        _dqdx.push_back((_calos[itcal].dQdx())[itHit]);
        if (correct_dQdx)
          _dqdx.at(_dqdx.size()-1) *= DqdxCorrection();
        _dedx.push_back((_calos[itcal].dEdx())[itHit]);
        _resrange.push_back((_calos[itcal].ResidualRange())[itHit]);
        _hitx.push_back(TrackPos.X());
    	  _hity.push_back(TrackPos.Y());
    	  _hitz.push_back(TrackPos.Z());
        _track_pitch.push_back((_calos[itcal].TrkPitchVec())[itHit]);
        if (correct_dQdx)
          _track_pitch.at(_track_pitch.size()-1) /= DqdxCorrection();
        // Get Hit Peak Time
        //const size_t & hitIndex = (_calos[itcal].TpIndices())[itHit];
        //const auto & thisHit = allHits[hitIndex];
        //hitPeakTime[planeNumb][itHit] = thisHit.PeakTime();
        _hitIndex.push_back((_calos[itcal].TpIndices())[itHit]);
        _hitPeakTime.push_back(INV_DBL); // For now in MCC11
        // Apply Lifetime corrections
        const geo::Point_t HitPoint(TrackPos.X(), TrackPos.Y(), TrackPos.Z());
        geo::TPCID const & tpcid = geom->FindTPCAtPosition(HitPoint);
        if (!tpcid.isValid) {
          std::cout << "CalorimetryHelper.cxx: " << "tpc not valid"<< std::endl;
          _drift_time.push_back(INV_DBL);
          _corr_factors.push_back(INV_DBL);
          continue;
        }
        int CryoID = geom->FindCryostatAtPosition(HitPoint);
        double Ticks = detprop->ConvertXToTicks(TrackPos.X(), planeNumb, tpcid.TPC, CryoID);
        _drift_time.push_back((Ticks - detprop->TriggerOffset()) * detprop->SamplingRate()*1e-3);
        _corr_factors.push_back(LifeTimeCorr(Ticks, 0));
      }
    }
    OrderResRange();
  }

  // Check if the calorimetry is valid
  bool CalorimetryHelper::IsValid() {

    _isValid = false;
    // Check that the size is not zero.
    if (_calos.size() == 0) _isValid = false;
    // If it's zero the loop is ignored
    for (size_t itcal = 0; itcal < _calos.size(); itcal++) {
      if (!(_calos[itcal].PlaneID().isValid)) continue;
      int planeNumb = _calos[itcal].PlaneID().Plane;
      if (_plane == planeNumb) {
        _isValid = true;
        break;
      }
    }
    return _isValid;

  }

  // Order the residual range with respect to the track direction
  void CalorimetryHelper::OrderResRange() {
    std::vector<double> res_vect;
      for (size_t i = 0; i < _resrange.size();i++) {
        res_vect.push_back(_resrange[i]);
      }
      _resrange_ord.resize(_resrange.size());
      size_t size = _resrange.size();
      double max = *max_element(res_vect.begin(),res_vect.end());
      for (size_t i = 0; i < _resrange.size();i++) {
        if (_hity[size-1] < _hity[0])  {
          if (_resrange[size-1] < _resrange[0])  {
            _resrange_ord[i] = _resrange[i];
          }
          else {
            _resrange_ord[i] = max-_resrange[i];
          }
        }
        else {
          if (_resrange[size-1] < _resrange[0])  {
            _resrange_ord[i] = max-_resrange[i];
          }
          else {
            _resrange_ord[i] = _resrange[i];
          }
        }
      }
  }

  // Get hit numb
  int CalorimetryHelper::GetHitNumb() {
    return _trackHitNumb;
  }
  // Get dqdx
  const std::vector<double> CalorimetryHelper::GetdQdx() {
    return _dqdx;
  }
  // Get dEdx
  const std::vector<double> CalorimetryHelper::GetdEdx() {
    return _dedx;
  }
  // Get Residual range
  const std::vector<double> CalorimetryHelper::GetResRange() {
    return _resrange;
  }
  // Get Ordered residual range
  const std::vector<double> CalorimetryHelper::GetResRangeOrdered() {
    return _resrange_ord;
  }
  // Get HitX
  const std::vector<double> CalorimetryHelper::GetHitX() {
    // std::cout << "hitx" << std::endl;
    // for (size_t i = 0; i < _hitx.size(); i++) {
    //   std::cout << _hitx[i] << std::endl;
    // }
    return _hitx;
  }
  // Get HitY
  const std::vector<double> CalorimetryHelper::GetHitY() {
    // std::cout << "hity" << std::endl;
    // for (size_t i = 0; i < _hity.size(); i++) {
    //   std::cout << _hity[i] << std::endl;
    // }
    return _hity;
  }
  // Get HitZ
  const std::vector<double> CalorimetryHelper::GetHitZ() {
    // std::cout << "hitz" << std::endl;
    // for (size_t i = 0; i < _hitz.size(); i++) {
    //   std::cout << _hitz[i] << std::endl;
    // }
    return _hitz;
  }
  // Get HitPeakTime
  const std::vector<double> CalorimetryHelper::GetHitPeakTime() {
    return _hitPeakTime;
  }
  // Get lifetime correction factors
  const std::vector<double> CalorimetryHelper::GetCorrFactor() {
    return _corr_factors;
  }
  // Get drift times
  const std::vector<double> CalorimetryHelper::GetDriftTime() {
    return _drift_time;
  }
  // Get track pitches
  const std::vector<double> CalorimetryHelper::GetTrackPitch() {
    return _track_pitch;
  }
  // Get hit indeces.
  const std::vector<size_t> CalorimetryHelper::GetHitIndex() {
    return _hitIndex;
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
  void CalorimetryHelper::FillHisto_dQdxVsRR(TH2D *h_dQdxVsRR) {
    const std::vector<double> &dQdx = GetdQdx();
    const std::vector<double> &resRangeOrd = GetResRangeOrdered();
    for (size_t it = 0; it < dQdx.size(); it++) {
      h_dQdxVsRR->Fill(resRangeOrd[it],dQdx[it]);
    }
  }

  // FIll 2D histo of dQdx vs residual range for hits in a given plane, in a track pitch interval
  void CalorimetryHelper::FillHisto_dQdxVsRR(TH2D *h_dQdxVsRR, const double &tp_min, const double &tp_max) {
    const std::vector<double> &dQdx = GetdQdx();
    const std::vector<double> &resRangeOrd = GetResRangeOrdered();
    const std::vector<double> &trackPitch = GetTrackPitch();
    for (size_t it = 0; it < dQdx.size(); it++) {
      if (trackPitch[it]<tp_min || trackPitch[it]>tp_max) continue;
      h_dQdxVsRR->Fill(resRangeOrd[it],dQdx[it]);
    }
  }

  // FIll 2D histo of dQdx vs residual range for hits in a given plane. Correct by MC lifetime
  void CalorimetryHelper::FillHisto_dQdxVsRR_LTCorr(TH2D *h_dQdxVsRR) {
    const std::vector<double> &dQdx = GetdQdx();
    const std::vector<double> &resRangeOrd = GetResRangeOrdered();
    const std::vector<double> &corrFactor = GetCorrFactor();
    for (size_t it = 0; it < dQdx.size(); it++)
      h_dQdxVsRR->Fill(resRangeOrd[it],dQdx[it]*corrFactor[it]);
  }

  // Same as above but with track pitch cut
  void CalorimetryHelper::FillHisto_dQdxVsRR_LTCorr(TH2D *h_dQdxVsRR, const double &tp_min, const double &tp_max) {
    const std::vector<double> &dQdx = GetdQdx();
    const std::vector<double> &resRangeOrd = GetResRangeOrdered();
    const std::vector<double> &trackPitch = GetTrackPitch();
    const std::vector<double> &corrFactor = GetCorrFactor();
    for (size_t it = 0; it < dQdx.size(); it++) {
      if (trackPitch[it]<tp_min || trackPitch[it]>tp_max) continue;
      h_dQdxVsRR->Fill(resRangeOrd[it],dQdx[it]*corrFactor[it]);
    }
  }

  // Fill 2D histo for dQdx/dEdx with lifetime correction. dEdx taken from MC.
  void CalorimetryHelper::FillHisto_dQdEVsRR_LTCorr_MC(TH2D *h_dQdEVsRR, const double &tp_min, const double &tp_max) {
    const std::vector<double> &dQdx = GetdQdx();
    const std::vector<double> &resRangeOrd = GetResRangeOrdered();
    const std::vector<double> &trackPitch = GetTrackPitch();
    const std::vector<double> &corrFactor = GetCorrFactor();
    for (size_t it = 0; it < dQdx.size(); it++) {
      if (trackPitch[it]<tp_min || trackPitch[it]>tp_max) continue;
      double dEdx = truedEdxHelper.GetMCdEdx(resRangeOrd[it]);
      h_dQdEVsRR->Fill(resRangeOrd[it],dQdx[it]*corrFactor[it]/dEdx);
    }
  }

  // Fill 2D histo for dQdx/dEdx with lifetime correction. dEdx taken from LandauVav.
  void CalorimetryHelper::FillHisto_dQdEVsRR_LTCorr_LV(TH2D *h_dQdEVsRR, const double &tp_min, const double &tp_max) {
    const std::vector<double> &dQdx = GetdQdx();
    const std::vector<double> &resRangeOrd = GetResRangeOrdered();
    const std::vector<double> &trackPitch = GetTrackPitch();
    const std::vector<double> &corrFactor = GetCorrFactor();
    for (size_t it = 0; it < dQdx.size(); it++) {
      if (trackPitch[it]<tp_min || trackPitch[it]>tp_max) continue;
      double rex = resRangeOrd[it];
      double dEdx = truedEdxHelper.LandauVav(rex,(tp_min+tp_max)/2);
      h_dQdEVsRR->Fill(resRangeOrd[it],dQdx[it]*corrFactor[it]/dEdx);
    }
  }

  // Set the parameters from the FHICL file
  void CalorimetryHelper::reconfigure(fhicl::ParameterSet const &p) {
    fTrackerTag = p.get<std::string>("TrackerTag");
    fCalorimetryTag = p.get<std::string>("CalorimetryTag");
    fPFParticleTag = p.get<std::string>("PFParticleTag");
  }

  // Return factor to correct dQdx values from bug in calorimetry module.
  double CalorimetryHelper::DqdxCorrection() {
    return geom->WirePitch(0)/geom->WirePitch(2);
  }

  // Reset
  void CalorimetryHelper::Reset() {
    _isValid = false;
    _isCalorimetrySet = false;
    _isData = false;
    _calos.clear();
    _trackHitNumb = INV_INT;
  	_dqdx.clear();
  	_dedx.clear();
  	_resrange.clear();
    _resrange_ord.clear();
    _hitx.clear();
  	_hity.clear();
  	_hitz.clear();
    _hitPeakTime.clear();
    _corr_factors.clear();
    _drift_time.clear();
    _track_pitch.clear();
    _hitIndex.clear();
  }

} // end of namespace stoppingcosmicmuonselection

#endif
