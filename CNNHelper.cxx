/***
  Class containing useful functions for geometry.

*/
#ifndef CNN_HELPER_CXX
#define CNN_HELPER_CXX

#include "CNNHelper.h"

namespace stoppingcosmicmuonselection {

  CNNHelper::CNNHelper() {

  }

  CNNHelper::~CNNHelper() {

  }

  // Get the score for a given hit.
  float CNNHelper::GetHitMichelScore(const anab::MVAReader<recob::Hit,4> &hitResults, const art::Ptr<recob::Hit> &hit)  {

    // Get CNN output for this hit.
    std::array<float,4> cnn_out = hitResults.getOutput(hit);
    // Get index of Michel output.
    size_t michelIndex = hitResults.getIndex("michel");
    // Get score.
    float score = cnn_out[michelIndex];

    return score;

  }

  // Get the number of michel hits according to a threshold.
  size_t CNNHelper::GetNumbMichelHits(const anab::MVAReader<recob::Hit,4> &hitResults, const artPtrHitVec &hits, float threshold) {

    size_t counter = 0;

    // Loop over hits.
    for (size_t i = 0; i < hits.size(); i++) {

      float score = GetHitMichelScore(hitResults, hits[i]);

      if (score > threshold) counter++;
    }

    return counter;

  }

  // Get the vector of scores.
  std::vector<double> CNNHelper::GetScoreVector(const anab::MVAReader<recob::Hit,4> &hitResults, const artPtrHitVec &hits) {

    std::vector<double> scores;
    scores.clear();

    for (size_t i = 0; i < hits.size(); i++) {

      // Get hit score.
      double score = GetHitMichelScore(hitResults, hits[i]);

      // Store score.
      scores.push_back(score);
    }

    return scores;

  }

  // Fill the 2D image of hits in the plane according to the score from the CNN.
  void CNNHelper::FillHitScoreImage(TProfile2D *image,
                                    const anab::MVAReader<recob::Hit,4> &hitResults,
                                    const artPtrHitVec &hits) {
    image->Reset();

    for (const art::Ptr<recob::Hit> &hitp : hits) {
      if (!hitp->WireID().isValid) continue;
      double hitPeakTime = hitp->PeakTime();
      unsigned int wireID = geoHelper.GetWireNumb(hitp);
      double score = GetHitMichelScore(hitResults,hitp);
      image->Fill(wireID,hitPeakTime,score);
    }
  }

  // Fill 1D histogram with the score for a given vector.
  void CNNHelper::FillScoreDistribution(TH1D *h, const anab::MVAReader<recob::Hit,4> &hitResults, const artPtrHitVec &hits) {

    h->Reset();

    if (hits.size() == 0) return;

    for (const art::Ptr<recob::Hit> &hitp : hits) {
      double score = GetHitMichelScore(hitResults, hitp);
      h->Fill(score);
    }
    return;
  }

} // end of namespace stoppingcosmicmuonselection

#endif
