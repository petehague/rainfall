#ifndef NUMBER_H
#define NUMBER_H

#include <cinttypes>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>

static uint16_t raninit_counter = 1;

class rng {
  std::mt19937 generator;
public:
  rng() {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    seed *= (unsigned)(raninit_counter++);
    generator.seed(seed);
  }

  double flat() {
    return (double)generator() / (double)UINT32_MAX;
  }

  double gaussian() {
    return sin(2.*M_PI*flat())*sqrt(-2*(log(generator())-log((double)UINT32_MAX)));;
  }
};

class unity {
  uint32_t i;
  bool overflow;
  rng rangen;

public:
  unity(double start);
  unity(uint32_t start);
  unity(uint32_t start, bool of);
  unity();

  ~unity();

  double val() const;
  uint32_t raw() const;
  bool OF() const;

  unity operator=(double x);
  unity operator=(uint32_t x);
  unity operator+(const unity &other);
  unity operator+(double x);
  unity operator-(const unity &other);
};



double meanvalue(std::vector <unity> numlist);
double stdev(std::vector <unity> numlist);

#endif
