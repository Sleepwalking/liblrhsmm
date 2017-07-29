/*
  liblrhsmm
  ===
  Copyright (c) 2016-2017 Kanru Hua. All rights reserved.

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

{
  int nseg = seg -> nseg;
  FP_TYPE* p_ = output_lg;
  FP_TYPE* a_;
# ifdef VITERBI
  if(forward != NULL) {
    *forward = calloc(nt * nseg, sizeof(FP_TYPE));
    a_ = *forward;
  } else
    a_ = calloc(nt * nseg, sizeof(FP_TYPE));
# else
  a_ = calloc(nt * nseg, sizeof(FP_TYPE));
# endif

# ifdef VITERBI
  int* s_ = calloc(nt * nseg, sizeof(int));
  int* reseg = calloc(nseg, sizeof(int));
  int* states = calloc(nt, sizeof(int));
  #define s(t, i) s_[(t) * nseg + (i)]
# endif

  for(int t = 0; t < nt; t ++)
    for(int i = 0; i < nseg; i ++)
      a(t, i) = NEGINF;

# ifdef VITERBI
  for(int i = 0; i < nseg; i ++)
    s(0, i) = 0;
# endif
  a(0, 0) = p(0, 0);

  int prune_range = lrh_inference_stprune_full_slope * nseg;
  for(int t = 1; t < nt; t ++) {
    for(int j = max(0, t * nseg / nt - prune_range);
            j < min(nseg, t * nseg / nt + prune_range); j ++) {
      FP_TYPE pstay = 1.0 - 1.0 / max(1.01, model -> durations[seg -> durstate[j]] -> mean);
      
#     ifdef VITERBI
      int maxfrom = j;
#     endif
      FP_TYPE maxlg = NEGINF;
      FP_TYPE p;

      // j -> j
      p = a(t - 1, j) + log(pstay) * lrh_daem_temperature;
#     ifdef VITERBI
      if(p > maxlg) {
        maxlg = p;
        maxfrom = j;
      }
#     else
      maxlg = lrh_lse(maxlg, p);
#     endif

      // j - 1 -> j
      if(j > 0) {
        FP_TYPE ptrans = 1.0 / max(1.01, model -> durations[seg -> durstate[j - 1]] -> mean);
        p = a(t - 1, j - 1) + log(ptrans) * lrh_daem_temperature;
#       ifdef VITERBI
        if(p > maxlg) {
          maxlg = p;
          maxfrom = j - 1;
        }
#       else
        maxlg = lrh_lse(maxlg, p);
#       endif
      }

      a(t, j) = maxlg + p(t, j);
#     ifdef VITERBI
      s(t, j) = maxfrom;
#     endif
    }
  }

  /*
  for(int t = 0; t < nt; t ++) {
    FP_TYPE maxt = -999999999;
    int pruned = lrh_inference_stprune_full_slope * nseg;
    for(int j = max(1, t * nseg / nt - pruned);
            j < min(nseg, t * nseg / nt + pruned); j ++)
      maxt = max(a(t, j), maxt);
    printf("%d %f %d\n", t, maxt, nseg);
  }
*/

# ifdef VITERBI
  int t = nt - 1;
  states[t] = nseg - 1;
  while(t > 0) {
    states[t - 1] = s(t, states[t]);
    t --;
  }

/*
  for(int t = 0; t < nt; t ++)
    printf(" %d", states[t]);
  puts("\n");
*/
  reseg[0] = 0;
  for(int t = 1; t < nt; t ++)
    if(states[t] != states[t - 1])
      reseg[states[t]] = t;

  // Shift
  for(int i = 0; i < nseg - 1; i ++)
    reseg[i] = reseg[i + 1];
  reseg[nseg - 1] = nt;
# endif

# ifdef VITERBI
  if(forward == NULL) free(a_);
  free(s_);
  free(states);
  return reseg;
# else
  return a_;
# endif
}
