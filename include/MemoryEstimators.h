#ifndef _MEMORY_ESTIMATORS_H
#define _MEMORY_ESTIMATORS_H

#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <list>
#include <numeric>
#include <cmath>

template <typename T, std::size_t N>
using VectorOfArrays = std::vector<std::array<T, N>>;

template <typename T, std::size_t N1, std::size_t N2>
using ArrayOfArrays = std::array<std::array<T, N1>, N2>;

template <typename T, std::size_t N>
using ListOfArrays = std::list<std::array<T, N>>;

template <typename T>
using VectorOfVectors = std::vector<std::vector<T>>;

template <typename T, std::size_t N>
T get_total_score(VectorOfArrays<T, N> segments, uint64_t buckets) {
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

template <typename T, std::size_t N1, std::size_t N2>
T get_total_score(ArrayOfArrays<T, N1, N2> segments, uint64_t buckets) {
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
size_t count_buckets(VectorOfArrays<T, N> segments) {
  auto bucket_count = 0; 
  for (auto it = segments.begin(); it != segments.end(); it++) {
    bucket_count += it->size();
  }
  return bucket_count;
}

template <typename T>
size_t count_buckets(VectorOfVectors<T> segments) {
  auto bucket_count = 0; 
  for (auto it = segments.begin(); it != segments.end(); it++) {
    bucket_count += it->size();
  }
  return bucket_count;
}

template <typename T, std::size_t N>
VectorOfVectors<T> keep_if_big_enough(VectorOfArrays<T, N> segments, T minimum_score) {
  VectorOfVectors<T> result;
  for (auto it = segments.begin(); it != segments.end(); it++) {
    std::vector<T> filtered;
    auto end = std::remove_if(it->begin(), it->end(), 
      [minimum_score](T score) { return score < minimum_score; });
    for (auto inner = it->begin(); inner != end; inner++) {
      filtered.push_back(*inner);
    }
    result.push_back(filtered);
  }
  return result;
}

template <typename T>
VectorOfVectors<T> keep_if_big_enough(VectorOfVectors<T> segments, T minimum_score) {
  VectorOfVectors<T> result;
  for (auto it = segments.begin(); it != segments.end(); it++) {
    std::vector<T> filtered;
    auto end = std::remove_if(it->begin(), it->end(), 
      [minimum_score](T score) { return score < minimum_score; });
    for (auto inner = it->begin(); inner != end; inner++) {
      filtered.push_back(*inner);
    }
    result.push_back(filtered);
  }
  return result;
}


template <typename T>
struct FragmentedResultIterator {
  using elem_type = T;
  typename std::vector<T>::const_iterator current;
  typename std::vector<T>::const_iterator end;
};

template <typename U>
std::vector<U> build_result_iterators(const VectorOfVectors<typename U::elem_type> &vector,
  size_t kFragmentSize) {
  typename std::vector<U> result;
  for (auto it = vector.begin(); it != vector.end(); it++) {
    for (auto score_it = it->cbegin(); score_it != it->cend(); score_it++) {
      unsigned long elements_seen = score_it - it->cbegin();
      if (elements_seen % kFragmentSize == 0) {
        unsigned long elementsLeft = it->cend() - score_it;
        auto iteratorSize = std::min(elementsLeft, kFragmentSize);
        result.emplace_back(U{score_it, score_it + iteratorSize});
      }
    }
  }
  return result;
}

inline uint8_t get_index_of_leftmost_one(uint64_t input) {
  if (input == 0) {
    return 65;
  }
  uint8_t index = 1;
  while ((input & (1UL << 63)) == 0) {
    input <<= 1;
    index++;
  }
  return index;
}

inline uint64_t get_upper(uint64_t size) {
  auto upper = 1UL << 32;
  while (upper < size) {
    upper <<= 1;
  }
  return upper;
}

inline double get_correction_factor(uint64_t size) {
  return 1.0 / 30;
}

template <typename T>
inline double calculate_indicator_function(std::vector<T> registers) {
  double sum = std::accumulate(registers.begin(), registers.end(), 0.0,
    [](double previous, T current) {
      return previous + pow(2, -1 * current);
    });
  return 1.0 / sum;
}

#endif // _MEMORY_ESTIMATORS_H