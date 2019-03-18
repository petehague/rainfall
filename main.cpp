/*
  The all new RainfallMCMC
*/


#include <iostream>

#include "include/rainfall.hpp"

using namespace std;

struct config {
  double targetRate = 0.25;
  uint32_t lookback = 100;
} settings;

extern double rainfall_init();
extern double log_likelihood(vector <atom*> p);
extern void per_model(model* m);
extern void rainfall_clear(vector <vector <atom*>> processQueue, vector <double *> results);
extern double rainfall_finish();

rng randomNumber;

double getrate(vector <bool> n) {
  double rate = 0;
  for (auto i=0;i<n.size();i++) {
    if (n[i]) { rate += 1; }
  }
  return rate/(double)n.size();
}

void resizeStep(vector <model*> chain, double *stepsize) {
  static vector <bool> acceptance;
  static uint32_t lbcounter = 0;
  static double rate;

  if (acceptance.size()<settings.lookback) {
    acceptance.push_back(chain.back()->lastaccept());
  } else {
    acceptance[lbcounter++] = chain.back()->lastaccept();
    if (lbcounter==acceptance.size()) { lbcounter = 0; }
    rate = getrate(acceptance);
    for (auto i=0;i<3;i++) {
      stepsize[i] *= (settings.targetRate-rate)/settings.targetRate;
      if (stepsize[i]<1e-6) { stepsize[i]=1e-6; }
      //cout << stepsize[i] << " ";
    }
  }
}

/*void resizeStep(vector <model*> chain, double *stepsize) {
  double flux = chain.back()->getAtom(0)->getParameter(2);
  for (auto i=0;i<2;i++) {
    stepsize[i] = 0.00017/(flux);
  }
}*/

int main() {
  vector <model*> chain;
  uint32_t maxmodels = 100000;
  bool verbose = true;
  model *currentModel;
  double stepsize[3];
  double lambda = 0;

  for (auto i=0;i<3;i++) { stepsize[i] = 1.; }

  rainfall_init();

  if (verbose) { cout << "Initiating Chain..." << endl; }
  chain.push_back(new model(3,1));
  chain[0]->compute(log_likelihood);

  for (auto i=0;i<maxmodels;i++) {
    if (verbose && i % (maxmodels/10) == 0) {
      cout << i << ": ";
      for (auto j=0;j<3;j++) {
        cout << " " << chain.back()->getAtom(0)->getParameter(j);
      }
      cout << " :" << chain.back()->llikelihood() << endl;
    }
    currentModel = chain.back()->step(log_likelihood, stepsize, lambda);
    currentModel->compute(log_likelihood);
    chain.push_back(currentModel);

    //resizeStep(chain, stepsize);

    if (i>maxmodels/2) { per_model(chain.back()); }
    else { lambda += 2./(double)maxmodels; }
  }

  for (auto i=0;i<chain.size();i++) {
    delete chain[i];
  }

  rainfall_finish();
}
