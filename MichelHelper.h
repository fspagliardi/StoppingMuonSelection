/***
  Class containing useful functions for Michel electron selection.

*/
#ifndef MICHEL_HELPER_H
#define MICHEL_HELPER_H

#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "TVector3.h"
#include "TMath.h"

#include "DataTypes.h"

namespace stoppingcosmicmuonselection {

  class MichelHelper {

  public:
    MichelHelper();
    ~MichelHelper();

    float GetHitMichelScore(const anab::MVAReader<recob::Hit,4> &hitResults, const art::Ptr<recob::Hit> &hit);

    size_t GetNumbMichelHits(const anab::MVAReader<recob::Hit,4> &hitResults, const artPtrHitVec &hits, float threshold);

    std::vector<double> GetScoreVector(const anab::MVAReader<recob::Hit,4> &hitResults, const artPtrHitVec &hits); 

  private:


  };
}

#endif
