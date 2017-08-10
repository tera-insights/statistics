#ifndef _MEMORY_ESTIMATORS_H
#define _MEMORY_ESTIMATORS_H

#include <vector>
#include <array>
#include <algorithm>
#include <iostream>

template <typename T, std::size_t N>
using VectorOfArrays = std::vector<std::array<T, N>>;

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
  typename std::vector<T>::const_iterator begin;
  const std::size_t number_of_elements;
};

template <typename T>
const std::vector<const FragmentedResultIterator<T>*> build_result_iterators(VectorOfVectors<T> vector,
  size_t kFragmentSize) {
  typename std::vector<const FragmentedResultIterator<T>*> result;
  for (auto it = vector.begin(); it != vector.end(); it++) {
    //std::cout << "it->size() = " << it->size() << std::endl;
    auto score_it = it->cbegin();
    for (;score_it != it->cend(); score_it++) {
      unsigned long elements_seen = score_it - it->cbegin();
      if (elements_seen % kFragmentSize == 0) {
        unsigned long elementsLeft = it->cend() - score_it;
        auto iteratorSize = std::min(elementsLeft, kFragmentSize);
        // std::cout << "iterator size = " << iteratorSize << std::endl;
        std::cout << "score_it = " << std::addressof(*score_it) << std::endl;
        std::cout << "*score_it = " << *score_it << std::endl;
        result.push_back(new FragmentedResultIterator<T>{score_it, iteratorSize});
        auto last_it = result.end() - 1;
        //std::cout << "most recently pushed fri = " << *((**last_it).begin) << std::endl;
      }
    }
  }
  return result;
}

#endif // _MEMORY_ESTIMATORS_H