#ifndef TOOLS_CXX
#define TOOLS_CXX

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include "TVector3.h"
#include "TMath.h"

#include "Tools.h"

namespace stoppingcosmicmuonselection {

  template<typename T>
  std::vector<std::vector<T>> get_neighbors(const std::vector<T> &object,
                                            const size_t &numbNeighbors) {

    if ((2*numbNeighbors+1)>object.size()) {
      std::cout << "Number of neighbors " << numbNeighbors
                << " is too big. Trying with "
                << numbNeighbors-1 << " neighbours." << std::endl;
      return get_neighbors(object, numbNeighbors-1);
    }
    std::vector<std::vector<T>> data;
    size_t objectSize = object.size();
    size_t m = numbNeighbors;

    for (size_t i = 0; i < objectSize; i++) {
      std::vector<T> inner;

      if (i < m) {
        for (size_t j =0; j <= 2*i; j++) {
          inner.push_back(object.at(j));
        }
      }
      else if (i > objectSize-m-1) {
        for (size_t j = 2*i-objectSize+1; j<objectSize;j++) {
          inner.push_back(object.at(j));
        }
      }
      else {
        for (size_t j = i-m; j <= i+m; j++) {
          inner.push_back(object.at(j));
        }
      }

      data.push_back(inner);

    }

    return data;
  }

  // Get the median without the max and min value.
  double get_smooth_trunc_median(std::vector<double> data) {
    if (data.size() > 2) {
      auto it_max = std::max_element(data.begin(), data.end());
      data.erase(it_max);

      auto it_min = std::min_element(data.begin(), data.end());
      data.erase(it_min);
    }

    double median = INV_DBL;

    size_t size = data.size();
    std::sort(data.begin(), data.end());
    if (size % 2 == 0) {
      median = (data[size/2 -1] + data[size/2]) / 2.;
    }
    else {
      median = data[size/2];
    }

    return median;
  }

  // Get the mean.
  double mean(const std::vector<double>& data)  {
    if (data.size() == 0)
      std::cout << "No data to calculate the mean." << std::endl;
    double result = 0;

    for (const auto & el : data)
      result += el;

    return (result / ((double)data.size()));
  }

  // Get the covariance.
  double cov (const std::vector<double>& data1,
              const std::vector<double>& data2) {
    if (data1.size()==0 || data2.size()==0)
      std::cout << "No data to calculate the covariance." << std::endl;
    if (data1.size() != data2.size())
      std::cout << "Data incompatible to calculate the covariance." << std::endl;

    double result = 0;
    auto mean1 = mean(data1);
    auto mean2 = mean(data2);

    for (size_t i = 0; i < data1.size(); i++) {
      result += (data1[i] - mean1)*(data2[i] - mean2);
    }

    result = result / ((double)data1.size());
    return result;
  }

  // Get the standard deviation.
  double stdev(const std::vector<double>& data) {
    if (data.size() == 0)
      std::cout << "No data to calculate the st. deviation." << std::endl;

    double result = 0;
    auto average = mean(data);

    for (auto const &el : data) {
      result += (el - average)*(el - average);
    }
    result = TMath::Sqrt(result / ((double)(data.size()-1)));
    return result;
  }

  // Print content of a vector.
  void printVec(const std::vector<double> &data) {
    for (size_t i = 0; i < data.size(); i++) {
      if (i == 0)
        std::cout << data[i];
      else {
        std::cout << ", " << data[i];
      }
    }
    std::cout << std::endl;
  }

}

#endif
