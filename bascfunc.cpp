#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "include/rainfall.hpp"

using namespace std;

double *dmap;
double *dbeam;
uint32_t naxis1;
uint32_t beamnaxis1;
uint32_t cra, cdec;
double noiseLevel;

inline uint32_t coords(uint32_t ra, uint32_t dec, uint32_t naxis) {
  return naxis*dec + ra;
}

double getNoise() {
  double mean, stdev;

  mean = 0;
  for (auto index=0;index<naxis1*naxis1;index++) {
    mean+=dmap[index];
  }
  mean /= (naxis1*naxis1);

  stdev = 0;
  for (auto index=0;index<naxis1*naxis1;index++) {
    stdev+=(dmap[index]-mean)*(dmap[index]-mean);
  }
  stdev /= (naxis1*naxis1)-1;

  return sqrt(stdev);
}

double log_likelihood(atom *params, uint8_t natoms) {
  double result = 0;

  for (auto i=0;i<natoms;i++) {
    uint32_t rabin_i = params[i].getParameter(0)*naxis1;
    uint32_t decbin_i = params[i].getParameter(1)*naxis1;
    double flux_i = params[i].getParameter(2);
    flux_i = noiseLevel*(flux_i/(1.-flux_i));
    result -= 2. * flux_i * dmap[coords(rabin_i, decbin_i,naxis1)];
    for (auto j=0;j<natoms;j++) {
      uint32_t dra = cra - (params[j].getParameter(0)*naxis1-rabin_i);
      uint32_t ddec = cdec - (params[j].getParameter(1)*naxis1-decbin_i);
      double flux_j = params[j].getParameter(2);
      flux_j = noiseLevel*(flux_j/(1.-flux_j));
      result += flux_i * flux_j * dbeam[coords(dra, ddec,beamnaxis1)];
    }
  }

  return -0.5*result/(noiseLevel*noiseLevel);
}

class bascfunc : public procQueue {
  vector <unity> x,y,F;
public:
  bascfunc() {
    fstream inputfile;

    naxis1 = 1024;
    beamnaxis1 = 1024;
    cra = 512;
    cdec = 512;

    dmap = new double[naxis1*naxis1];
    dbeam = new double[naxis1*naxis1*4];

    inputfile.open("ex_image.txt");
    for (auto ra=0;ra<naxis1;ra++) {
      for (auto dec=0;dec<naxis1;dec++) {
        inputfile >> dmap[coords(ra,dec,naxis1)];
      }
    }
    inputfile.close();

    inputfile.open("ex_psf.txt");
    for (auto ra=0;ra<beamnaxis1;ra++) {
      for (auto dec=0;dec<beamnaxis1;dec++) {
        inputfile >> dbeam[coords(ra,dec,beamnaxis1)];
      }
    }
    inputfile.close();

    noiseLevel = getNoise();

    cout << "Noise Level: " << noiseLevel << endl;
  }

  void clear() {
    indices.push_back(bufferSize);
    for (auto i=0;i<indices.size()-1;i++) {
      uint32_t start = indices[i];
      uint32_t length = indices[i+1]-start;
      *result[i] = log_likelihood(&buffer[start], length);
    }
  }

  void tick(model *m) {
    x.push_back(m->getAtom(0)->getUnity(0));
    y.push_back(m->getAtom(0)->getUnity(1));
    F.push_back(m->getAtom(0)->getUnity(2));
  }

  ~bascfunc() {
    double meanf, upper, lower, stdf;

    cout << "x: " << meanvalue(x)*naxis1 << " +- " << stdev(x)*naxis1 << endl;
    cout << "y: " << meanvalue(y)*naxis1 << " +- " << stdev(y)*naxis1 << endl;

    meanf = meanvalue(F);
    meanf = noiseLevel*(meanf/(1.-meanf));
    upper = meanvalue(F)+stdev(F);
    lower = upper - 2*stdev(F);
    upper = noiseLevel*(upper/(1.-upper));
    lower = noiseLevel*(lower/(1.-lower));
    stdf = 0.5*(upper-lower);
    cout << "F: " << meanf << " +- " << stdf << endl;

    delete[] dmap;
    delete[] dbeam;
  }
};

bascfunc *evaluator = new bascfunc;
