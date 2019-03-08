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

double rainfall_init() {
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

double log_likelihood(vector <atom*> params) {
  uint8_t natoms = params.size();
  double result = 0;

  for (auto i=0;i<natoms;i++) {
    uint32_t rabin_i = params[i]->getParameter(0)*naxis1;
    uint32_t decbin_i = params[i]->getParameter(1)*naxis1;
    double flux_i = params[i]->getParameter(2);
    result += 2. * flux_i * dmap[coords(rabin_i, decbin_i,naxis1)];
    for (auto j=0;j<natoms;j++) {
      uint32_t dra = cra - (params[j]->getParameter(0)*naxis1-rabin_i);
      uint32_t ddec = cdec - (params[j]->getParameter(1)*naxis1-decbin_i);
      double flux_j = params[j]->getParameter(2);
      result -= flux_i * flux_j * dbeam[coords(dra, ddec,beamnaxis1)];
    }
  }

  return -0.5*result/(noiseLevel*noiseLevel);
}

double rainfall_finish() {

}
