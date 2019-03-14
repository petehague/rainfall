/*
  The all new RainfallMCMC
*/


#include <iostream>

#include "include/rainfall.hpp"

using namespace std;

extern double rainfall_init();
extern double log_likelihood(vector <atom*> p);
extern void per_model(model* m);
extern double rainfall_finish();

rng randomNumber;

double getrate(vector <bool> n) {
  double rate = 0;
  for (auto i=0;i<n.size();i++) {
    if (n[i]) { rate += 1; }
  }
  return rate/(double)n.size();
}

int main() {
  vector <model*> chain;
  uint32_t maxmodels = 10000;
  bool verbose = true;
  model *currentModel;
  double stepsize[3];
  vector <bool> acceptance;
  uint32_t lookback = 100;
  uint32_t lbcounter = 0;
  double rate, targetRate = 0.5;

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
    currentModel = chain.back()->step(log_likelihood, stepsize);
    currentModel->compute(log_likelihood);
    chain.push_back(currentModel);
    if (acceptance.size()<lookback) {
      acceptance.push_back(chain.back()->lastaccept());
    } else {
      acceptance[lbcounter++] = chain.back()->lastaccept();
      if (lbcounter==acceptance.size()) { lbcounter = 0; }
      rate = getrate(acceptance);
      for (auto i=0;i<3;i++) {
        stepsize[i] *= (targetRate-rate)/targetRate;
        if (stepsize[i]<1e-6) { stepsize[i]=1e-6; }
        //cout << stepsize[i] << " ";
      }
    }

    per_model(chain.back());
  }

  for (auto i=0;i<chain.size();i++) {
    delete chain[i];
  }

  rainfall_finish();
}
