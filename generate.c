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

#include "generate.h"
#include "math-funcs.h"

lrh_observ* lrh_generate(lrh_model* src, lrh_seg* seq) {
  for(int i = 0; i < seq -> nseg; i ++) {
    lrh_duration* dr = src -> durations[seq -> durstate[i]];
    int idur = round(fmax(1, lrh_gaussian_random(dr -> mean, dr -> var)));
    if(i == 0) seq -> time[i] = idur;
    else seq -> time[i] = seq -> time[i - 1] + idur;
  }
  int nt = seq -> time[seq -> nseg - 1];
  int* dims = calloc(src -> nstream, sizeof(int));
  for(int i = 0; i < src -> nstream; i ++)
    dims[i] = src -> streams[i] -> gmms[0] -> ndim;
  
  lrh_observ* ret = lrh_create_observ(src -> nstream, nt, dims);
  for(int i = 0; i < seq -> nseg; i ++) {
    for(int t = i == 0 ? 0 : seq -> time[i - 1]; t < seq -> time[i]; t ++) {
      for(int l = 0; l < src -> nstream; l ++) {
        lrh_gmm* pdf = src -> streams[l] -> gmms[seq -> outstate[l][i]];
        int mixture = lrh_random_choose(pdf -> weight, pdf -> nmix);
        for(int k = 0; k < dims[l]; k ++)
          lrh_obm(ret, t, k, l) = lrh_gaussian_random(lrh_gmmu(pdf, mixture, k),
            lrh_gmmv(pdf, mixture, k));
      }
    }
  }
  
  free(dims);
  return ret;
}

