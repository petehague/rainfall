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

extern procQueue *evaluator;

rng randomNumber;

int main() {
  vector <model*> chain;
  uint32_t maxmodels = 10;
  bool verbose = true;
  model *currentModel;
  double lambda = 0;

  if (verbose) { cout << "Initiating Chain..." << endl; }
  chain.push_back(new model(3,1));
  chain[0]->compute(evaluator);

  for (auto i=0;i<maxmodels;i++) {
    if (verbose && i % (maxmodels/10) == 0) {
      cout << i << ": ";
      for (auto j=0;j<3;j++) {
        cout << " " << chain.back()->getAtom(0)->getParameter(j);
      }
      cout << " :" << chain.back()->llikelihood() << endl;
    }
    currentModel = chain.back()->step(evaluator, lambda);
    currentModel->compute(evaluator);
    chain.push_back(currentModel);

    if (i>=maxmodels/2) { evaluator->tick(currentModel); }
    else { lambda += 2./(double)maxmodels; }
  }

  delete evaluator;

  for (auto i=0;i<chain.size();i++) {
    delete chain[i];
  }
}
