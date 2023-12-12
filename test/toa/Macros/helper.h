#ifndef Macros_helper_h
#define Macros_helper_h

// C headers
#include <numeric>


std::vector<size_t> decrease_sorted_indices(const std::vector<double>& v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indices based on comparing values in v (decreasing order)
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  return idx;
};
  

//time-interval based on that ~210ps wide and with the highest number of hits
//extension valid in high PU of taking smallest interval with (order of)68% of hits
std::vector<float> fixSizeHighestDensity(const std::vector<double>& time, std::vector<double> weight, const unsigned int& minNhits=3, const float& deltaT=0.210, const float& timeWidthBy=0.5)
{
  if (time.size() < minNhits)
    return std::vector<float>({-99., -1., -99});

  if (weight.empty())
    weight.resize(time.size(), 1.);

  std::vector<float> t(time.size(), 0.);
  std::vector<float> w(time.size(), 0.);
  std::vector<size_t> sortedIndex = decrease_sorted_indices(time);
  for (std::size_t i = 0; i < sortedIndex.size(); ++i) {
    t[i] = time[sortedIndex[i]];
    w[i] = weight[sortedIndex[i]];
  }

  int max_elements = 0;
  int start_el = 0;
  int end_el = 0;
  float timeW = 0.f;
  float tolerance = 0.05f;

  for (auto start = t.begin(); start != t.end(); ++start) {
    const auto startRef = *start;
    int c = count_if(start, t.end(), [&](float el) { return el - startRef <= deltaT + tolerance; });
    if (c > max_elements) {
      max_elements = c;
      auto last_el = find_if_not(start, t.end(), [&](float el) { return el - startRef <= deltaT + tolerance; });
      auto valTostartDiff = *(--last_el) - startRef;
      if (std::abs(deltaT - valTostartDiff) < tolerance) {
        tolerance = std::abs(deltaT - valTostartDiff);
      }
      start_el = distance(t.begin(), start);
      end_el = distance(t.begin(), last_el);
      timeW = valTostartDiff;
    }
  }

  // further adjust time width around the chosen one based on the hits density
  // proved to improve the resolution: get as many hits as possible provided they are close in time
  float HalfTimeDiff = timeW * timeWidthBy;
  float sum = 0.;
  float num = 0;
  int totSize = t.size();
  
  for (int ij = 0; ij <= start_el; ++ij) {
    if (t[ij] > (t[start_el] - HalfTimeDiff)) {
      for (int kl = ij; kl < totSize; ++kl) {
        if (t[kl] < (t[end_el] + HalfTimeDiff)) {
          sum += t[kl] * w[kl];
          num += w[kl];
        } else
          break;
      }
      break;
    }
  }

  if (num == 0)
    return std::vector<float>({-99., -1., -99});
  
  return std::vector<float>({sum / num, float(1. / sqrt(num)), float(totSize), t[start_el] - HalfTimeDiff, t[end_el] + HalfTimeDiff});
};


#endif // #ifndef Macros_helper_h
