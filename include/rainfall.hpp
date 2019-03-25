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
  atom();
  atom(uint8_t n);
  ~atom();
  atom(const atom& other);
  atom& operator=(const atom& other);
  uint8_t getSize() const { return size; }
  double getParameter(uint8_t n) const { return params[n].val(); }
  unity getUnity(uint8_t n) const { return params[n]; }
  void setParameter(uint8_t n, double value) { params[n] = value; }
  atom step();
};

class model;

class procQueue {
protected:
  atom *buffer; //Flat array of atom data
  uint32_t bufferSize;
  vector <double*> result; //Where the log likelihoods will go
  vector <uint32_t> indices; //Indices in buffer of the starts of each model

public:
  procQueue();
  virtual ~procQueue();
  void push (vector <atom> model, double *returnLoc);
  virtual void clear();
  virtual void tick(model *m) = 0;
  void pop();
};

class model {
  vector <atom> atoms;
  double loglike;
  bool accepted;

public:
  model();
  model(uint8_t nparams, uint8_t natoms);
  ~model();
  atom *getAtom(uint8_t atom_n);
  void compute(procQueue* evaluator);
  double llikelihood() const { return loglike; }
  bool lastaccept() const { return accepted; }
  model* step(procQueue* evaluator, double lambda);
};

extern rng randomNumber;

#endif
