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
  const artPtrHitVec HitHelper::GetArtPtrToHitVect(const art::FindManyP<recob::Hit> &fmht,
                                                   const size_t &trackIndex) {
    return fmht.at(trackIndex);
  }

  // Get hit list on a given plane.
  artPtrHitVec HitHelper::GetHitsOnAPlane(const size_t &planeNumb,
                                          const artPtrHitVec &allHits) {
    artPtrHitVec hitsOnPlane;
    for (auto const &hitp : allHits) {
      if (!hitp->WireID().isValid) continue;
      if (hitp->WireID().Plane != planeNumb) continue;
      hitsOnPlane.push_back(hitp);
    }
    return hitsOnPlane;
  }

  // Get XYZ position for that hit
  TVector3 HitHelper::GetHitXYZ(const art::Ptr<recob::Hit> &hitp,
                                art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                const std::vector<art::Ptr<recob::Track>> &tracklist,
                                const size_t &trackIndex) {

    TVector3 hitLoc(INV_DBL,INV_DBL,INV_DBL);
    if (!fmthm.isValid()) return hitLoc;

    const artPtrHitVec vhit = fmthm.at(trackIndex);
    const std::vector<const recob::TrackHitMeta*> vmeta = fmthm.data(trackIndex);
    // iterate on meta data
    for (size_t ii=0;ii<vhit.size();++ii) {
      if (vhit[ii].key() != hitp.key())
        continue;
      if (vmeta[ii]->Index() == std::numeric_limits<int>::max()) {
        continue;
      }
      if (vmeta[ii]->Index()>=tracklist[trackIndex]->NumberTrajectoryPoints()){
        throw cet::exception("HitHelper.cxx") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[trackIndex]->NumberTrajectoryPoints()<<" for track index "<<trackIndex<<". Something is wrong";
      }
      if (!tracklist[trackIndex]->HasValidPoint(vmeta[ii]->Index())){
        continue;
      }
      hitLoc = tracklist[trackIndex]->LocationAtPoint<TVector3>(vmeta[ii]->Index());
    } // iteration on metadata

    return hitLoc;
  }

  // Get index of the closest hit to a given point on the given track.
  const size_t HitHelper::GetIndexClosestHitToPoint(const TVector3 &point,
                                                    const artPtrHitVec &hits,
                                                    art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                                    const std::vector<art::Ptr<recob::Track>> &tracklist,
                                                    const size_t &trackIndex) {
    std::vector<double> distanceVec;
    for (auto const &hitp : hits) {
      const TVector3 &hitLoc = GetHitXYZ(hitp,fmthm,tracklist,trackIndex);
      distanceVec.push_back((hitLoc - point).Mag());
    }
    auto it_min_element = std::min_element(distanceVec.begin(),distanceVec.end());
    size_t hitIndex = it_min_element - distanceVec.begin();
    std::cout << "Closest hit to point: " << std::endl << "\tWireID: "
              << geoHelper.GetWireNumb(hits[hitIndex])
              << "\tTime: " << hits[hitIndex]->PeakTime() << std::endl;
    return hitIndex;
  }

  // Get the closest hit to a given point on the given track.
  art::Ptr<recob::Hit> HitHelper::GetClosestHitToPoint(const TVector3 &point,
                                                       const artPtrHitVec &hits,
                                                       art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                                       const std::vector<art::Ptr<recob::Track>> &tracklist,
                                                       const size_t &trackIndex) {
    size_t hitIndex = GetIndexClosestHitToPoint(point, hits, fmthm, tracklist, trackIndex);
    return hits.at(hitIndex);
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

  // Get subvector of michel-like hits.
  artPtrHitVec HitHelper::GetMichelLikeHits(const artPtrHitVec &hits,
                                            const TVector3 &recoEndPoint,
                                            art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                            const std::vector<art::Ptr<recob::Track>> &tracklist,
                                            const size_t &trackIndex) {

    std::vector<art::Ptr<recob::Hit>> result;
    result.clear();

    for (const art::Ptr<recob::Hit> &hitp : hits) {
      if (!IsHitMichelLike(hitp,recoEndPoint,fmthm,tracklist,trackIndex))
        continue;
      result.push_back(hitp);
    }

    if (result.size() == 0)
      std::cout << "HitHelper.cxx: " << "HitHelper::GetMichelLikeHits " << "is returning an empty vector." << std::endl;
    return result;
  }

  // Get subvector of muon-like hits.
  artPtrHitVec HitHelper::GetMuonLikeHits(const artPtrHitVec &hits,
                                          const TVector3 &recoEndPoint,
                                          art::FindManyP<recob::Hit,recob::TrackHitMeta> &fmthm,
                                          const std::vector<art::Ptr<recob::Track>> &tracklist,
                                          const size_t &trackIndex) {

    std::vector<art::Ptr<recob::Hit>> result;
    result.clear();

    for (const art::Ptr<recob::Hit> &hitp : hits) {
      if (IsHitMichelLike(hitp,recoEndPoint,fmthm,tracklist,trackIndex))
        continue;
      result.push_back(hitp);
    }

    if (result.size() == 0)
      std::cout << "HitHelper.cxx: " << "HitHelper::GetMuonlLikeHits " << "is returning an empty vector." << std::endl;
    return result;
  }

  // Fill the TGraph2D for the images.
  void HitHelper::FillTrackGraph2D(TGraph2D *graph,
                                   const artPtrHitVec &trackHits,
                                   const TVector3 &recoEndPoint,
                                   const size_t &planeNumber) {
    graph->Set(0);
    const geo::GeometryCore *geom = lar::providerFrom<geo::Geometry>();
    const geo::Point_t EndPoint(recoEndPoint.X(), recoEndPoint.Y(), recoEndPoint.Z());
    geo::TPCID const & tpcid = geom->FindTPCAtPosition(EndPoint);
    if (!tpcid.isValid) {
      std::cout << "Track End Point is in invalid TPC. Return empty Graph." << std::endl;
      return;
    }
    for (size_t i = 0; i < trackHits.size(); i++) {
      const art::Ptr<recob::Hit> &hitp = trackHits[i];
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
      size_t wireOffset = geoHelper.GetWireOffset(hit_tpcid, planeNumber);
      graph->SetPoint(i,wireID+wireOffset,hitPeakTime,electron_perc);
    }
    return;
  }

  // Get a TProfile2D filled with hit peak times and wire number
  void HitHelper::FillTrackHitPicture(TProfile2D* image,
                                     const artPtrHitVec &trackHits,
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
      size_t wireOffset = geoHelper.GetWireOffset(hit_tpcid, planeNumber);
      image->Fill(wireID+wireOffset,hitPeakTime,electron_perc);
    }
    return;
  }

  // Initialise the image for a series of hit for a given plane
  void HitHelper::InitHitImageHisto(TProfile2D *&image, const size_t &planeNumber, const std::string &name) {
    size_t nWires = geoHelper.GetNumberWiresOneSide(planeNumber);
    int nTicks = detprop->NumberTimeSamples();
    image = new TProfile2D(name.c_str(),name.c_str(),nWires,0,nWires,nTicks,0,nTicks);
  }

  // Set the parameters from the FHICL file
  void HitHelper::reconfigure(fhicl::ParameterSet const &p) {
    _electronEnergyFractionToCallMichelHits = p.get<double>("electronEnergyFractionToCallMichelHits", 0.7);
    _maxDistanceToCallMichelHits = p.get<double>("maxDistanceToCallMichelHits", 15);
  }

} // end of namespace stoppingcosmicmuonselection

#endif
