#include "../include/MemoryEstimators.h"
#include <array>
#include <iostream>
#include <cassert>

void test_get_total_score_on_ints() {
  std::vector<std::array<int, 3>> entries;
  entries.push_back(std::array<int, 3>{1, 2, 3});
  entries.push_back(std::array<int, 3>{4, 5, 6});
  assert(get_total_score(entries, 6) == 21);
  assert(get_total_score(entries, 5) == 15);
}

void test_count_buckets() {
  std::vector<std::array<int, 3>> entries;
  assert(count_buckets(entries) == 0);
  std::array<int, 3> entryOne = {1, 2, 3};
  entries.push_back(entryOne);
  assert(count_buckets(entries) == 3);
  entries.push_back(entryOne);
  assert(count_buckets(entries) == 6);
}

void test_keep_if_big_enough_100_percent_multiplier() {
  std::vector<std::array<int, 3>> entries;
  entries.push_back(std::array<int, 3>{1, 2, 3});
  entries.push_back(std::array<int, 3>{4, 5, 6});
  auto bucket_count = count_buckets(entries);
  auto filtered = keep_if_big_enough(entries, 7);
  bucket_count = count_buckets(filtered);
  assert(bucket_count == 0);
}

void test_keep_if_big_enough_keep_one_element_then_filter_again() {
  std::vector<std::array<int, 3>> entries;
  entries.push_back(std::array<int, 3>{1, 2, 3});
  entries.push_back(std::array<int, 3>{4, 5, 6});
  auto bucket_count = count_buckets(entries);
  auto filtered = keep_if_big_enough(entries, 6);
  bucket_count = count_buckets(filtered);
  assert(bucket_count == 1);
  assert(filtered.at(1).size() == 1);
  assert(filtered.at(1).at(0) == 6);
  filtered = keep_if_big_enough(filtered, 6);
  bucket_count = count_buckets(filtered);
  assert(bucket_count == 1);
  assert(filtered.at(1).size() == 1);
  assert(filtered.at(1).at(0) == 6);
}

void test_build_result_iterators() {
  std::vector<std::vector<int>> vector;
  vector.push_back(std::vector<int>{1, 2, 3});
  vector.push_back(std::vector<int>{4, 5, 6});
  auto result = build_result_iterators(vector, 2);
  int elements_seen = 0;
  for (size_t i = 0; i < result.size() - 1; i++) {
    for (size_t j = 0; j < result.at(i)->number_of_elements; j++) {
      auto iterator = result.at(i)->begin + j;
      size_t outer_index = elements_seen / 3;
      size_t inner_index = elements_seen - (3 * outer_index);
      //std::cout << "vector.at(outer_index).at(inner_index) = " << vector.at
      //(outer_index).at(inner_index) << std::endl;
      std::cout << "iterator = " << std::addressof(*iterator) << std::endl;
      std::cout << "*iterator = " << *iterator << std::endl;
      assert(vector.at(outer_index).at(inner_index) == *iterator);
      elements_seen++;
    }
  }
}

int main() {
  test_get_total_score_on_ints();
  test_count_buckets();
  test_keep_if_big_enough_100_percent_multiplier();
  test_keep_if_big_enough_keep_one_element_then_filter_again();
  test_build_result_iterators();
  std::cout << "All good!" << std::endl;
  return 0;
}