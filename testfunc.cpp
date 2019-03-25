#include <vector>
#include <iostream>

#include "include/rainfall.hpp"

using namespace std;

class testfunc : public procQueue {
public:
  void clear() {
    for (auto i=0;i<indices.size();i++) {
      *result[i] = 0.5;
    }
  }

  void tick(model *) {

  }
};

testfunc *evaluator = new testfunc;
