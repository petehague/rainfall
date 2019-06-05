#include "include/rainfall.hpp"

model::model() { }

model::model(uint8_t nparams, uint8_t natoms) {
  for (auto i=0;i<natoms;i++) {
    atoms.push_back(atom(nparams));
    stepsize.push_back(0.5);
  }
}

model::~model() {
  atoms.clear();
}

atom* model::getAtom(uint8_t atom_n) {
  return &atoms[atom_n];
}

void model::compute(procQueue* evaluator) {
  double result;
  evaluator->push(atoms, &result);
  evaluator->clear();
  evaluator->pop();
  loglike = result;
}

double* model::prepare_step(procQueue *evaluator) {
  double *gradpoints = new double[atoms[0].getSize()*2];

  for (auto i=0;i<atoms[0].getSize();i++) {
    vector <atom> proposal;
    double dx = randomNumber.gaussian()*stepsize[i];
    for (auto k=-1;k<2;k+=2) {
      for (auto j=0;j<atoms[0].getSize();j++) {
        if (i==j) {
          proposal.push_back(atoms[0].getParameter(j)+dx*(double)k);
        } else {
          proposal.push_back(atoms[0].getParameter(j));
        }
      }
    }
    evaluator->push(proposal,&gradpoints[i]);
  }

  return gradpoints;
}

model* model::resolve_step(procQueue *evaluator, double lambda, double gradpoints) {
  uint16_t resultIndex = 0;

  for (auto i=0;i<atoms[0].getSize();i++) {

    resultIndex++;
  }

  //delete[] gradpoints;
}

model* model::step(procQueue *evaluator, double lambda) {
  model* destination = new model(atoms[0].getSize(),atoms.size());
  double newlike;

  for (auto i=0;i<atoms.size();i++) {
    for (auto j=0;j<atoms[0].getSize();j++) {
      double dx = randomNumber.gaussian();
      destination->atoms[i].setParameter(j, atoms[i].getParameter(j)+dx);
    }
  }

  compute(evaluator);
  destination->compute(evaluator);
  newlike=destination->llikelihood();

  accepted = true;

  /*
   Rejection critera; make the new model the same as the previous one
   Use this expression as division is safer than subtraction in floating point
  */
  if (exp(newlike*lambda)/exp(loglike*lambda)<randomNumber.flat()) {
  //if (newlike<loglike) {
    for (auto i=0;i<atoms.size();i++) {
      for (auto j=0;j<atoms[0].getSize();j++) {
        double dx = randomNumber.gaussian();
        destination->atoms[i].setParameter(j, atoms[i].getParameter(j));
      }
    }
    accepted = false;
  }
  return destination;
}
