#include "include/rainfall.hpp"

atom::atom() {
  size = 0;
}

atom::atom(uint8_t n) {
  params = new unity[n];
  size = n;
}

atom::~atom() {
  if (size>0 && params!=NULL) {  delete[] params; }
}

atom::atom(const atom &other) {
  params = new unity[other.size];
  size = other.size;
}

atom& atom::operator=(const atom& other) {
  if (size>0 && params!=NULL) { delete[] params; }
  size = other.size;
  params = new unity[size];
  for (auto i=0;i<size;i++) {
    params[i] = other.params[i];
  }
  return *this;
}

atom atom::step() {
  return atom(size);
}
