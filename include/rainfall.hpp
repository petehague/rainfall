#ifndef RAINFALL_H
#define RAINFALL_H

#include <cinttypes>
#include <cinttypes>
#include <vector>

#include "number.hpp"

class atom {
  unity *params;
  uint8_t size;

public:
  atom() {
    size = 0;
  }

  atom(uint8_t n) {
    params = new unity[n];
    size = n;
  }

  ~atom() {
    if (size>0 && params!=NULL) {  delete[] params; }
  }

  uint8_t getSize() const { return size; }
  double getParameter(uint8_t n) const { return params[n].val(); }
  unity getUnity(uint8_t n) const { return params[n]; }
  void setParameter(uint8_t n, double value) {
    params[n] = value;
  }


  atom& operator=(atom other) {
    if (size>0 && params!=NULL) { delete[] params; }
    size = other.size;
    params = new unity[size];
    for (auto i=0;i<size;i++) {
      params[i] = other.params[i];
    }
    return *this;
  }

  atom step() {
    return atom(size);
  }
};

#endif
