/***
  Class containing useful functions for geometry.

*/
#ifndef MICHEL_HELPER_CXX
#define MICHEL_HELPER_CXX

#include "MichelHelper.h"

namespace stoppingcosmicmuonselection {

  MichelHelper::MichelHelper() {

  }

  MichelHelper::~MichelHelper() {

  }

  float MichelHelper::GetHitMichelScore(const anab::MVAReader<recob::Hit,4> &hitResults, const art::Ptr<recob::Hit> &hit)  {
    
    // Get CNN output for this hit.
    std::array<float,4> cnn_out = hitResults.getOutput(hit);
    // Get index of Michel output.
    size_t michelIndex = hitResults.getIndex("michel");
    // Get score.
    float score = cnn_out[michelIndex];

    return score;

  }

  size_t MichelHelper::GetNumbMichelHits(const anab::MVAReader<recob::Hit,4> &hitResults, const artPtrHitVec &hits, float threshold) {
    
    size_t counter = 0;

    // Loop over hits.
    for (size_t i = 0; i < hits.size(); i++) {
      
      float score = GetHitMichelScore(hitResults, hits[i]);

      if (score > threshold) counter++;
    }

    return counter;

  }

} // end of namespace stoppingcosmicmuonselection

#endif
