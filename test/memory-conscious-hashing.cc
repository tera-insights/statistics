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
  std::cout << "got total score on ints" << std::endl;
}

void test_get_number_of_buckets() {
  std::vector<std::array<int, 3>> entries;
  assert(get_number_of_buckets(entries) == 0);
  std::array<int, 3> entryOne = {1, 2, 3};
  entries.push_back(entryOne);
  assert(get_number_of_buckets(entries) == 3);
  entries.push_back(entryOne);
  assert(get_number_of_buckets(entries) == 6);
  std::cout << "got total number of buckets" << std::endl;
}

int main() {
  test_get_total_score_on_ints();
  test_get_number_of_buckets();
  std::cout << "Tests passed" << std::endl;
  return 0;
}