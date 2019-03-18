#ifndef RAINFALL_H
#define RAINFALL_H

#include <cinttypes>
#include <cinttypes>
#include <vector>

#include "number.hpp"

using namespace std;

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

extern rng randomNumber;

class model {
  vector <atom*> atoms;
  double loglike;
  bool accepted;

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

  bool lastaccept() const { return accepted; }

  model* step(double (*likefunc)(vector <atom*>), double *stepsize, double lambda) {
    model* destination = new model(atoms[0]->getSize(),atoms.size());
    double newlike;

    for (auto i=0;i<atoms.size();i++) {
      for (auto j=0;j<atoms[0]->getSize();j++) {
        double dx = randomNumber.gaussian()*stepsize[j];
        destination->atoms[i]->setParameter(j, atoms[i]->getParameter(j)+dx);
      }
    }

    compute(likefunc);
    destination->compute(likefunc);
    newlike=destination->llikelihood();

    accepted = true;

    /*
     Rejection critera; make the new model the same as the previous one
     Use this expression as division is safer than subtraction in floating point
    */
    if (exp(newlike*lambda)/exp(loglike*lambda)<randomNumber.flat()) {
    //if (newlike<loglike) {
      for (auto i=0;i<atoms.size();i++) {
        for (auto j=0;j<atoms[0]->getSize();j++) {
          double dx = randomNumber.gaussian();
          destination->atoms[i]->setParameter(j, atoms[i]->getParameter(j));
        }
      }
      accepted = false;
    }
    return destination;
  }
};

#endif
