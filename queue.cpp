
#include "include/rainfall.hpp"
#include <iostream>

#define MAX_ATOMS 10000

procQueue::procQueue() {
  buffer = new atom[MAX_ATOMS];
}

procQueue::~procQueue() {
  delete[] buffer;
}

void procQueue::push (vector <atom> model, double *returnLoc) {
  if (bufferSize+model.size()>MAX_ATOMS-1) { return; }

  indices.push_back(bufferSize);
  for (auto i=bufferSize;i<bufferSize+model.size();i++) {
    buffer[i] = model[i-bufferSize]; //Note we want to copy the atom itself
  }
  bufferSize += model.size();

  result.push_back(returnLoc);
}

void procQueue::pop() {
  bufferSize = 0;
  result.clear();
  indices.clear();
}

void procQueue::clear() {
}

//void procQueue::tick(model *m) {
//}
