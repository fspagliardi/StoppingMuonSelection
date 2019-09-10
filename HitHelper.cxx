/***
  Class containing useful functions for hit.

*/
#ifndef HIT_HELPER_CXX
#define HIT_HELPER_CXX

#include "HitHelper.h"

namespace stoppingcosmicmuonselection {

  HitHelper::HitHelper() {

  }

  HitHelper::~HitHelper() {

  }

  // Get the track index of a track object
  size_t HitHelper::GetTrackIndex(const recob::Track &track,
                                  const std::vector<art::Ptr<recob::Track>> &tracklist) {
    if (tracklist.size() == 0) return INV_INT;
    for (size_t trackIter = 0; trackIter < tracklist.size(); trackIter++) {
      if (tracklist[trackIter]->ID() == track.ID())
        return trackIter;
    }
    return INV_INT;
  }

  // Get the vector of art::Ptr to hit for the given track
  const std::vector<art::Ptr<recob::Hit>> HitHelper::GetArtPtrToHitVect(const art::FindManyP<recob::Hit> &fmht,
                                                                        const size_t &trackIndex) {
    return fmht.at(trackIndex);
  }

  // Get XYZ position for that hit
  TVector3 HitHelper::GetHitXYZ(const art::Ptr<recob::Hit> &hitp,
                                art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                const std::vector<art::Ptr<recob::Track>> &tracklist,
                                const size_t &trackIndex) {

    TVector3 hitLoc(INV_DBL,INV_DBL,INV_DBL);
    if (!fmthm.isValid()) return hitLoc;

    const std::vector<art::Ptr<recob::Hit>> vhit = fmthm.at(trackIndex);
    const std::vector<const recob::TrackHitMeta*> vmeta = fmthm.data(trackIndex);
    // iterate on meta data
    for (size_t ii=0;ii<vhit.size();++ii) {
      if (vhit[ii].key() != hitp.key())
        continue;
      if (vmeta[ii]->Index() == std::numeric_limits<int>::max()) {
        continue;
      }
      if (vmeta[ii]->Index()>=tracklist[trackIndex]->NumberTrajectoryPoints()){
        throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[trackIndex]->NumberTrajectoryPoints()<<" for track index "<<trackIndex<<". Something is wrong";
      }
      if (!tracklist[trackIndex]->HasValidPoint(vmeta[ii]->Index())){
        continue;
      }
      hitLoc = tracklist[trackIndex]->LocationAtPoint<TVector3>(vmeta[ii]->Index());
    } // iteration on metadata

    return hitLoc;
  }

  // Check if a hit has high electron contribution at a certain distance from the end point
  bool HitHelper::IsHitMichelLike(const art::Ptr<recob::Hit> &hitp,
                                  const TVector3 &recoEndPoint,
                                  art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                  const std::vector<art::Ptr<recob::Track>> &tracklist,
                                  const size_t &trackIndex) {
    for (const sim::TrackIDE &tIDE : bt_serv->HitToTrackIDEs(*hitp)) {
      // check if there is a contribution from an electron.
      if (TMath::Abs(pi_serv->TrackIdToParticle_P(tIDE.trackID)->PdgCode())==11) {
        TVector3 hitLoc = GetHitXYZ(hitp,fmthm,tracklist,trackIndex);
        if (hitLoc == TVector3(INV_DBL,INV_DBL,INV_DBL)) continue;
        if (tIDE.energyFrac>_electronEnergyFractionToCallMichelHits
            && (hitLoc-recoEndPoint).Mag()<_maxDistanceToCallMichelHits)
          return true;
      }
    }
    return false;
  }

  // Get a TProfile2D filled with hit peak times and wire number
  void HitHelper::FillTrackHitPicture(TProfile2D* image,
                                     const std::vector<art::Ptr<recob::Hit>> &trackHits,
                                     const TVector3 &recoEndPoint,
                                     const size_t &planeNumber) {
    image->Reset();
    const geo::GeometryCore *geom = lar::providerFrom<geo::Geometry>();
    const geo::Point_t EndPoint(recoEndPoint.X(), recoEndPoint.Y(), recoEndPoint.Z());
    geo::TPCID const & tpcid = geom->FindTPCAtPosition(EndPoint);
    if (!tpcid.isValid) {
      std::cout << "Track End Point is in invalid TPC. Filling the image at (1,1)." << std::endl;
      // Dummy filling
      image->Fill(1,1,1);
      return;
    }
    size_t nWires = geoHelper.GetNumberWiresOneSide(planeNumber);
    for (const art::Ptr<recob::Hit> & hitp : trackHits) {
      // Get only hit in the collection plane
      if (!hitp->WireID().isValid) continue;
      if (hitp->WireID().Plane != planeNumber) continue;
      unsigned int hit_tpcid = hitp->WireID().TPC;
      double hitPeakTime = hitp->PeakTime();
      unsigned int wireID = hitp->WireID().Wire;
      double electron_perc = 0;
      for(const sim::TrackIDE& ide : bt_serv->HitToTrackIDEs(*hitp)) {
        if (TMath::Abs(pi_serv->TrackIdToParticle_P(ide.trackID)->PdgCode())==11) {//contribution from electron
          electron_perc += ide.energyFrac;
          //std::cout << "Electron perc: " << electron_perc << std::endl;
        }
      }
      if (hit_tpcid == geoHelper.tpcIndecesBL[0] || hit_tpcid == geoHelper.tpcIndecesBR[0])
        image->Fill(wireID,hitPeakTime,electron_perc);
      else if (hit_tpcid == geoHelper.tpcIndecesBL[1] || hit_tpcid == geoHelper.tpcIndecesBR[1])
        image->Fill((nWires/3+wireID),hitPeakTime,electron_perc);
      else if (hit_tpcid == geoHelper.tpcIndecesBL[2] || hit_tpcid == geoHelper.tpcIndecesBR[2])
        image->Fill((2*nWires/3+wireID),hitPeakTime,electron_perc);
    }
    return;
  }

  // Initialise the image for a series of hit for a given plane
  void HitHelper::InitHitImageHisto(TProfile2D *image, const size_t &planeNumber, const std::string &name) {
    size_t nWires = geoHelper.GetNumberWiresOneSide(planeNumber);
    int nTicks = detprop->NumberTimeSamples();
    image = new TProfile2D(name.c_str(),name.c_str(),nWires,0,nWires,nTicks,0,nTicks);
  }

  // Set the parameters from the FHICL file
  void HitHelper::reconfigure(fhicl::ParameterSet const &p) {
    fTrackerTag = p.get<std::string>("TrackerTag");
    fPFParticleTag = p.get<std::string>("PFParticleTag");
    _electronEnergyFractionToCallMichelHits = p.get<double>("electronEnergyFractionToCallMichelHits", 0.7);
    _maxDistanceToCallMichelHits = p.get<double>("maxDistanceToCallMichelHits", 15);
  }

} // end of namespace stoppingcosmicmuonselection

#endif
