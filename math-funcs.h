/*
  liblrhsmm
  ===
  Copyright (c) 2016 Kanru Hua. All rights reserved.

  This file is part of liblrhsmm.

  liblrhsmm is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  liblrhsmm is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with liblrhsmm.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef LRHSMM_MATH_FUNCS_H
#define LRHSMM_MATH_FUNCS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "external/fastapprox/fasttrig.h"
#include "external/fastapprox/fastlog.h"
#include "external/fastapprox/fastexp.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EPS (1e-7)
#define NEGINF (-1e15)

#define lrh_lse_table_size 7000
#define lrh_lse_table_step 0.001
#define lrh_exp_table_size 7000
#define lrh_exp_table_step 0.002
extern FP_TYPE lrh_lse_table[];
extern FP_TYPE lrh_exp_table[];

inline static FP_TYPE lrh_random() {
  return ((FP_TYPE)rand()) / RAND_MAX;
}

inline static FP_TYPE lrh_gaussian_random(FP_TYPE mean, FP_TYPE var) {
  FP_TYPE u1 = lrh_random();
  FP_TYPE u2 = lrh_random();
  return sqrt(-2.0 * log(u1) * var) * fastcos(2.0 * M_PI * u2) + mean;
}

inline static int lrh_random_choose(FP_TYPE* pmf, int n) {
  FP_TYPE x = lrh_random();
  FP_TYPE a = 0;
  for(int i = 0; i < n; i ++) {
    a += pmf[i];
    if(a > x) return i;
  }
  return n - 1;
}

// log-sum-exp function
inline static FP_TYPE lrh_lse(FP_TYPE a, FP_TYPE b) {
  FP_TYPE maxab = a > b ? a : b;
  FP_TYPE minab = a > b ? b : a;
  if(maxab - minab > 6.91) // if a and b differ by 3 orders of magnitude
    return maxab;
  return lrh_lse_table[(int)((maxab - minab) / lrh_lse_table_step)] + maxab;
  //return fastlog(1.0 + fastexp(minab - maxab)) + maxab;
}

inline static FP_TYPE lrh_negexp(FP_TYPE x) {
  if(x < -13.9) return 0.0;
  if(x >= 0) return 1.0;
  return lrh_exp_table[(int)((-x) / lrh_exp_table_step)];
  //return exp(x);
}

inline static FP_TYPE lrh_exp(FP_TYPE x) {
  if(x < -13.9) return 0.0;
  if(x >= 0) return fastexp(x);
  return lrh_exp_table[(int)((-x) / lrh_exp_table_step)];
  //return exp(x);
}

#endif

