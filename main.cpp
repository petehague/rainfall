/*
  The all new RainfallMCMC
*/


#include <iostream>

#include "include/rainfall.hpp"

using namespace std;

extern double rainfall_init();
extern double log_likelihood(vector <atom*> p);
extern double rainfall_finish();

rng randomNumber;

class model {
  vector <atom*> atoms;
  double loglike;

public:
  model() { }

  model(uint8_t nparams, uint8_t natoms) {
    for (auto i=0;i<natoms;i++) {
      atoms.push_back(new atom(nparams));
    }
  }

  ~model() {
    for (auto i=0;i<atoms.size();i++) {
      delete atoms[i];
    }
  }

  atom *getAtom(uint8_t atom_n) {
    return atoms[atom_n];
  }

  void compute(double (*likefunc)(vector <atom*>)) {
    loglike = likefunc(atoms);
  }

  double llikelihood() const { return loglike; }

  model* step(double (*likefunc)(vector <atom*>)) {
    model* destination = new model(atoms[0]->getSize(),atoms.size());
    double newlike;

    for (auto i=0;i<atoms.size();i++) {
      for (auto j=0;j<atoms[0]->getSize();j++) {
        double dx = randomNumber.gaussian();
        destination->atoms[i]->setParameter(j, atoms[i]->getParameter(j)+dx);
      }
    }

    compute(likefunc);
    destination->compute(likefunc);
    newlike=destination->llikelihood();

    if ((newlike-loglike)<log(randomNumber.flat())) {
      for (auto i=0;i<atoms.size();i++) {
        for (auto j=0;j<atoms[0]->getSize();j++) {
          double dx = randomNumber.gaussian();
          destination->atoms[i]->setParameter(j, atoms[i]->getParameter(j));
        }
      }
    }
    return destination;
  }
};

int main() {
  vector <model*> chain;
  uint32_t maxmodels = 1000;
  bool verbose = true;
  model *currentModel;

  rainfall_init();

  if (verbose) { cout << "Initiating Chain..." << endl; }
  chain.push_back(new model(3,1));
  chain[0]->compute(log_likelihood);

  for (auto i=0;i<maxmodels;i++) {
    if (verbose && i % 100 == 0) {
      cout << i << ": ";
      for (auto j=0;j<3;j++) {
        cout << " " << chain.back()->getAtom(0)->getParameter(j);
      }
      cout << " :" << chain.back()->llikelihood() << endl;
    }
    currentModel = chain.back()->step(log_likelihood);
    currentModel->compute(log_likelihood);
    chain.push_back(currentModel);
  }

  for (auto i=0;i<chain.size();i++) {
    delete chain[i];
  }
}
