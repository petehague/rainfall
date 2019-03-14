#include "include/number.hpp"

unity::unity(double start) {
  i = (uint32_t)(start/(double)UINT32_MAX);
  overflow = false;
}

unity::unity(uint32_t start) {
  i = start;
  overflow = false;
}

unity::unity(uint32_t start, bool of) {
  i = start;
  overflow = of;
}

unity::unity() {
  i = (uint32_t)(rangen.flat()*(double)UINT32_MAX);
}

unity::~unity() {
}

double unity::val() const {
  return (double)i/(double)UINT32_MAX;
}

uint32_t unity::raw() const {
  return i;
}

bool unity::OF() const {
  return overflow;
}

unity unity::operator+(const unity &other) {
  int64_t a = i;
  int64_t b = other.raw();

  if ((a+b)>(int64_t)UINT32_MAX) {
    return unity((uint32_t)UINT32_MAX, true);
  }
  return unity((uint32_t)(a+b));
}

unity unity::operator-(const unity &other) {
  int64_t a = i;
  int64_t b = other.raw();

  if ((a-b)<0) {
    return unity((uint32_t)0, true);
  }
  return unity((uint32_t)(a-b));
}

unity unity::operator=(uint32_t x) {
  i = x;
  return *this;
}

unity unity::operator=(double x) {
  i = (uint32_t)(x*(double)UINT32_MAX);
  return *this;
}

double meanvalue(std::vector <unity> numlist) {
  double result;

  result = 0;
  for (auto index=0;index<numlist.size();index++) {
    result+=numlist[index].val();
  }
  result /= numlist.size();

  return result;
}

double stdev(std::vector <unity> numlist) {
  double mean, result;

  mean = meanvalue(numlist);

  result = 0;
  for (auto index=0;index<numlist.size();index++) {
    result+=(numlist[index].val()-mean)*(numlist[index].val()-mean);
  }
  result /= numlist.size()-1;

  return sqrt(result);
}
