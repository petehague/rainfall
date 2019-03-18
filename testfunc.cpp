#include <vector>

#include "include/rainfall.hpp"

using namespace std;

double rainfall_init() {

}

double log_likelihood(vector <atom*> params) {
  return 0.5;
}

void rainfall_clear(vector <vector <atom*>> processQueue, vector <double *> results) {
  for (auto i=0;i<processQueue.size();i++) {
    results[i] = log_likelihood(processQueue[i]);
  }
}


void per_model(model *m) {

}

double rainfall_finish() {

}
