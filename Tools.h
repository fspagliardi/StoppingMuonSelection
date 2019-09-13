#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "TVector3.h"

#include "DataTypes.h"

namespace stoppingcosmicmuonselection {

  template<typename T>
  std::vector<std::vector<T>> get_windows(const std::vector<T>& the_thing,
                                        const size_t window_size);

  // Get the median without the max and min value.
  double get_smooth_trunc_median(std::vector<double> data);

  // Get the mean.
  double mean(const std::vector<double>& data);


  // Get the covariance.
  double cov (const std::vector<double>& data1,
              const std::vector<double>& data2);

  // Get the standard deviation.
  double stdev(const std::vector<double>& data);

}

#endif
