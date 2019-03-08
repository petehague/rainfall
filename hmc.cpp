/*

  Implementation of Hamiltonian Monte Carlo
  From http://www.mcmchandbook.net/HandbookChapter5.pdf

*/

#define MAXCHAIN 100000
#define RANSEED 1234

#include <iostream>
#include <cinttypes>
#include <random>
#include <cmath>

#include "include/pick.hpp"

using namespace std;

generator rangen(RANSEED);

double *dmap;
double *dbeam;
double *dmap_gradx;
double *dmap_grady;
double *dbeam_gradx;
double *dbeam_grady;
int32_t mapsize = 512;

double chain[MAXCHAIN];
double *chainptr[MAXCHAIN];

inline uint32_t coord(double x, double y, uint16_t s) {
  uint16_t ix = (uint16_t)x;
  uint16_t iy = (uint16_t)y;
  return ix+iy*s;
}

void scan(double *data, uint16_t size) {
  uint16_t block=16;

  for(auto y=0;y<size/block;y++) {
    for(auto x=0;x<size/block;x++) {
      double value = 0;
      for (auto v=y*block;v<(y+1)*block;v++) {
        for (auto u=x*block;u<(x+1)*block;u++) {
          value += data[u+v*size];
        }
      }
      cout << " ";
      if (value>0.1) { cout << "#"; continue; }
      if (value>0.01) { cout << "*"; continue; }
      if (value>0.001) { cout << "."; continue; }
      cout << " ";
    }
    cout << endl;
  }
}

void makemap(uint16_t mx, uint16_t my) {
  double dx,dy,r1sq,r2sq;
  uint32_t position;
  double sigmasq = 100;


  dmap = new double[mx*my];
  dbeam = new double[mx*my*4];
  dmap_gradx = new double[mx*my];
  dbeam_gradx = new double[mx*my*4];
  dmap_grady = new double[mx*my];
  dbeam_grady = new double[mx*my*4];
  
  //Dirty Map
  for (auto y=0;y<my;y++) {
    for (auto x=0;x<mx;x++) {
      dx = (double)(x - mapsize/2);
      dy = (double)(y - mapsize/2);
      r1sq = dx*dx + dy*dy;
      //dx = (double)(30 + x - mapsize/2);
      //dy = (double)(y - mapsize/2);
      //r2sq = dx*dx + dy*dy;
      dmap[coord(x,y,mx)] = exp(-r1sq/sigmasq);// + exp(-r2sq/sigmasq); 
    }
  }

  //Dirty map gradient
  for (auto y=0;y<my;y++) {
    for (auto x=0;x<mx;x++) {
      position = coord(x,y,mx);
      if (x==0 || y==0 || x==mx-1 || y==my-1) {
        dmap_gradx[position] = 0;
        dmap_grady[position] = 0;       
      } else {
        double M = 0.5*(-dmap[position-1]+dmap[position+1]);
        dmap_gradx[position] = 2. * M * x; 
        M = 0.5*(dmap[position-mx]+dmap[position+mx]);
        dmap_grady[position] = 2. * M * y;
      }
    }
  }

  //Dirty Beam  
  for (auto y=0;y<2*my;y++) {
    for (auto x=0;x<2*mx;x++) {
      dx = (double)(x - mapsize/2);
      dy = (double)(y - mapsize/2);
      r1sq = dx*dx + dy*dy;
      dbeam[coord(x,y,mx*2)] = exp(-r1sq/sigmasq); 
    }
  }

  //Dirty Beam gradient 
  for (auto y=0;y<2*my;y++) {
    for (auto x=0;x<2*mx;x++) {
      position = coord(x,y,2*mx);
      if (x==0 || y==0 || x==2*mx-1 || y==2*my-1) {
        dbeam_gradx[position] = 0;
        dbeam_grady[position] = 0;       
      } else {
        double M = 0.5*(-dbeam[position-1]+dbeam[position+1]);
        dbeam_gradx[position] = 2. * M * x; 
        M = 0.5*(dbeam[position-2*mx]+dbeam[position+2*mx]);
        dbeam_grady[position] = 2. * M * y;
      }
    }
  }
}

double loglike(double *model, uint16_t modelsize) {
  uint16_t natoms = modelsize/3;
  double ll = 0;  
  double Xn,Yn,Fn,Xm,Ym,Fm;

  for (uint16_t n=0;n<natoms;n++) {
     double *modelptr1 = model+n*3;
     Xn = modelptr1[0]*(double)mapsize;
     Yn = modelptr1[1]*(double)mapsize;
     Fn = modelptr1[2]*(1-modelptr1[2]);
     ll += 2.*Fn*dmap[coord(Xn, Yn, mapsize)];
     for (uint16_t m=0;m<natoms;m++) {
       double *modelptr2 = model+m*3;
       Xm = modelptr2[0]*(double)mapsize;
       Ym = modelptr2[1]*(double)mapsize;
       Fm = modelptr2[2]*(1-modelptr2[2]);
       uint32_t position = coord(mapsize + Xn - Xm, 
                                 mapsize + Yn - Ym,
                                 mapsize*2);
       ll -= Fn*Fm*dbeam[position];
     }
  }  

  return ll;
}

void llikegrad(double *model, uint16_t modelsize, double *grad) {
  uint16_t natoms = modelsize/3;
  double Xn,Yn,Fn,Xm,Ym,Fm;


  grad[0] = 0;
  grad[1] = 0;
  grad[2] = 0;

  // F component
  for (uint16_t n=0;n<natoms;n++) {
     double *modelptr1 = model+n*3;
     Xn = modelptr1[0]*(double)mapsize;
     Yn = modelptr1[1]*(double)mapsize;
     grad[0] += 2.*dmap[coord(Xn, Yn, mapsize)];
     for (uint16_t m=0;m<natoms;m++) {
       double *modelptr2 = model+m*3;
       Xm = modelptr2[0]*(double)mapsize;
       Ym = modelptr2[1]*(double)mapsize;
       Fm = modelptr2[2]*(1-modelptr2[2]);
       uint32_t position = coord(mapsize + Xn - Xm,
                                 mapsize + Yn - Ym,
                                 mapsize*2);
       grad[0] -= Fm*dbeam[position];
     }
  }  

  //xy component
  for (uint16_t n=0;n<natoms;n++) {
    double *modelptr1 = model+n*3;
    Xn = modelptr1[0]*(double)mapsize;
    Yn = modelptr1[1]*(double)mapsize;
    Fn = modelptr1[2]*(1-modelptr1[2]);
    grad[1] += 2.*Fn*dmap_gradx[coord(Xn, Yn, mapsize)];
    grad[2] += 2.*Fn*dmap_grady[coord(Xn, Yn, mapsize)];
    for (uint16_t m=0;m<natoms;m++) {
      double *modelptr2 = model+m*3;
      Xm = modelptr2[0]*(double)mapsize;
      Ym = modelptr2[1]*(double)mapsize;
      Fm = modelptr2[2]*(1-modelptr2[2]);
      uint32_t position = coord(mapsize + Xn - Xm, 
                                mapsize + Yn - Ym,
                                mapsize*2);
      grad[1] -= Fn*Fm*dbeam_gradx[position];
      grad[2] -= Fn*Fm*dbeam_grady[position]; 
    }
  }
}

inline void bounce(double *x, uint16_t n) {
  for (auto i=0;i<n;i++) {
    while (x[i]<0 || x[i]>1) {
      if (x[i]<0) { x[i]++; }
      if (x[i]>1) { x[i]--; }
    }
  }
}

double *hmc_step(double (U)(double *x, uint16_t n),
               void (grad_U)(double *x, uint16_t n, double *g),
               double epsilon,
               double L, 
               double *cur_q,
               uint16_t qcount) {

  double q[qcount];
  double p[qcount];
  double cur_p[qcount];
  double gradient[qcount];
  double current_U, current_K, proposed_U, proposed_K;

  rangen.getNormBlock(p, qcount);
  for (auto i=0;i<qcount;i++) { 
    q[i] = cur_q[i]; 
    cur_p[i] = p[i];
  }
 
  //Half step
  grad_U(q, qcount, gradient);
  for (auto i=0;i<qcount;i++) {
    p[i] -= epsilon*gradient[i]/2.;
  }

  //Alternating full steps
  for (auto step=0;step<L-1;step++) {
    for (auto i=0;i<qcount;i++) { q[i] += p[i]*epsilon; }
    bounce(q, qcount);
    grad_U(q, qcount, gradient);
    for (auto i=0;i<qcount;i++) { p[i] -= gradient[i]*epsilon; }
  }

  //Momentum half step at the end
  for (auto i=0;i<qcount;i++) { q[i] += p[i]*epsilon; }
  bounce(q, qcount);
  grad_U(q, qcount, gradient);
  for (auto i=0;i<qcount;i++) { 
    p[i] -= gradient[i]*0.5*epsilon; 
    p[i] = -p[i];
  }

  current_U = U(cur_q, qcount);
  current_K = 0;
  for (auto i=0;i<qcount;i++) { current_K += cur_p[i]*cur_p[i]*0.5; }
  proposed_U = U(q, qcount);
  proposed_K = 0;
  for (auto i=0;i<qcount;i++) { proposed_K += p[i]*p[i]*0.5; }

  if (rangen.getFlat()<exp(current_U-proposed_U+current_K-proposed_K)) {
    for (auto i=0;i<qcount;i++) { cur_q[i] = q[i]; }
  }
  return cur_q; 
}

int main() {
  uint16_t natoms;

  natoms = 1;
  makemap(mapsize, mapsize); 

  scan(dmap, mapsize);
  //scan(dbeam, mapsize*2);

  chainptr[0] = &chain[0];
  for (auto i=0;i<natoms;i++) {
    chain[i] = rangen.getFlat()*0.5+0.25;
  }

  for (auto i=0;i<1000;i++) {
    hmc_step(loglike, llikegrad, 0.1, 1.0, chainptr[i], natoms*3);
    chainptr[i+1] = chainptr[i]+3*natoms;
    for (auto j=0;j<natoms*3;j++) {
      cout << "  " << *(chainptr[i]+j);
    }
    cout << endl;
  }
}
