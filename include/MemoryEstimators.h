#ifndef _MEMORY_ESTIMATORS_H
#define _MEMORY_ESTIMATORS_H

#include <vector>
#include <array>

template <typename T, std::size_t N>
T get_total_score(std::vector<std::array<T, N>> segments, uint64_t buckets) {
  T total_score = 0;
  auto buckets_seen = 0;
  for (auto it = segments.begin(); it != segments.end(); it++) {
    auto segment_size = it->size();
    for (size_t i = 0; i < segment_size && buckets_seen < buckets; i++) {
      total_score += it->at(i);
      buckets_seen++;
    }
  }
  return total_score;
}

template <typename T, std::size_t N>
size_t get_number_of_buckets(std::vector<std::array<T, N>> segments) {
  auto bucket_count = 0; 
  for (auto it = segments.begin(); it != segments.end(); it++) {
    bucket_count += it->size();
  }
  return bucket_count;
}

#endif // _MEMORY_ESTIMATORS_H