#include "../include/MemoryEstimators.h"
#include <array>
#include <iostream>
#include <cassert>

void test_get_total_score_on_ints() {
  std::array<int, 3> entryOne = {1, 2, 3};
  std::array<int, 3> entryTwo = {4, 5, 6};
  std::vector<std::array<int, 3>> entries;
  entries.push_back(entryOne);
  entries.push_back(entryTwo);
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
  std::array<int, 3> entryOne = {1, 2, 3};
  std::array<int, 3> entryTwo = {4, 5, 6};
  std::vector<std::array<int, 3>> entries;
  entries.push_back(entryOne);
  entries.push_back(entryTwo);
  auto bucket_count = count_buckets(entries);
  auto filtered = keep_if_big_enough(entries, 7);
  bucket_count = count_buckets(filtered);
  assert(bucket_count == 0);
}

void test_keep_if_big_enough_keep_one_element_then_filter_again() {
  std::array<int, 3> entryOne = {1, 2, 3};
  std::array<int, 3> entryTwo = {4, 5, 6};
  std::vector<std::array<int, 3>> entries;
  entries.push_back(entryOne);
  entries.push_back(entryTwo);
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

int main() {
  test_get_total_score_on_ints();
  test_count_buckets();
  test_keep_if_big_enough_100_percent_multiplier();
  test_keep_if_big_enough_keep_one_element_then_filter_again();
  std::cout << "All good!" << std::endl;
  return 0;
}