#ifndef SV_H
#define SV_H
#include <cmath>
#include <iostream>
#include <utility>
using namespace std;
#include "error.hpp"
#include "Vector.hpp"
#include "full_mat_c.hpp"
#include "abstract_mat_c.hpp"
#include "mesh.hpp"
typedef double (*pfn) (double);


bool SV(int i, int itermax);

#endif // VOLFINI_H
