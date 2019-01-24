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
  FP_TYPE* a_begin = calloc(nseg * nt, sizeof(FP_TYPE));
  FP_TYPE* a_end   = calloc(nseg * nt, sizeof(FP_TYPE));

# ifdef VITERBI
  int* s_begin = calloc(nt * nseg, sizeof(int));
  int* s_end   = calloc(nt * nseg, sizeof(int));
  int* t_reseg = calloc(nseg, sizeof(int));
  int* s_reseg = calloc(nseg, sizeof(int));
  int  n_reseg = 0;
# define sb(t, i) s_begin[(t) * nseg + (i)]
# define se(t, i) s_end  [(t) * nseg + (i)]
# endif

# define ab(t, i) a_begin[(t) * nseg + (i)]
# define ae(t, i) a_end  [(t) * nseg + (i)]
  for(int t = 0; t < nt; t ++)
  for(int i = 0; i < nseg; i ++) {
    ab(t, i) = NEGINF;
    ae(t, i) = NEGINF;
  }
  for_tj_forward(0, 0, 1, 1)
    int jdur = seg -> time[j] - (j == 0 ? 0 : seg -> time[j - 1]);
    int maxdur, mindur = 1;
    lrh_duration* srcdr = model -> durations[seg -> durstate[j]];
    calculate_maxdur(srcdr, jdur);
    
#   ifdef VITERBI
    int maxfrom = max(0, t - jdur);
#   endif

    // compute initiation probability
    FP_TYPE pinit = NEGINF;
    int i = 0;
    do { // for each transition into the j-th state
      int srcst = j + seg -> djump_in[j][i];
      if(srcst >= 0 && t > 0) {
        FP_TYPE pj = ae(t - 1, srcst) +
          log(seg -> pjump_in[j][i]) * lrh_daem_temperature;
#       ifdef VITERBI
        if(pj > pinit) {
          pinit = pj;
          maxfrom = srcst;
        }
#       else
        pinit = lrh_lse(pinit, pj);
#       endif
      }
      i ++;
    } while(seg -> djump_in[j][i - 1] != -1);
    pinit += p(t, j);
    if(j == 0 && t > 0) pinit = NEGINF;
    else if(t == 0 && j > 0) pinit = NEGINF;
    else if(t == 0 && j == 0) pinit = p(t, j);
    
    ab(t, j) = pinit;
#   ifdef VITERBI
    sb(t, j) = maxfrom;
#   endif

    // compute termination probability
    FP_TYPE pterm = NEGINF;
    FP_TYPE pobserv = 0;
    for(int d = 1; d < mindur; d ++) {
      int t0 = t - d + 1;
      if(t0 < 0) break;
      pobserv += p(t0, j);
    }
    for(int d = mindur; d < maxdur; d ++) {
      int t0 = t - d + 1;
      if(t0 < 0) break;
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, d);
      FP_TYPE porigin = ab(t0, j);
#     ifdef VITERBI
      if(pobserv + pdur + porigin > pterm) {
        pterm = pobserv + pdur + porigin;
        maxfrom = t0;
      }
#     else
      pterm = lrh_lse(pterm, pobserv + pdur + porigin);
#     endif
      pobserv += p(t0, j);
    }
    ae(t, j) = pterm;
#   ifdef VITERBI
    se(t, j) = maxfrom;
#   endif

  end_for_tj()

# ifdef VITERBI
  // finalization
  s_reseg[0] = nseg - 1;
  t_reseg[0] = se(nt - 1, nseg - 1);
  n_reseg ++;

  // backtracking
  int i = s_reseg[0];
  int t = t_reseg[0];
  while(i != 0 && t > 0) {
    i = sb(t, i);
    t = se(t - 1, i);
    s_reseg = realloc(s_reseg, sizeof(int) * (n_reseg + 1));
    t_reseg = realloc(t_reseg, sizeof(int) * (n_reseg + 1));
    s_reseg[n_reseg] = i;
    t_reseg[n_reseg] = t;
    n_reseg ++;
  }
  
  // merge and reverse segmentation sequences
  int* reseg = calloc(n_reseg * 2 + 2, sizeof(int));
  for(int i = 0; i < n_reseg; i ++) {
    reseg[i * 2 + 0] = i == n_reseg - 1 ? nt : t_reseg[n_reseg - i - 2];
    reseg[i * 2 + 1] = s_reseg[n_reseg - i - 1];
  }
  reseg[n_reseg * 2 + 0] = -1;
  reseg[n_reseg * 2 + 1] = -1;
# endif

# undef ab
# undef ae
# ifdef VITERBI
  if(forward != NULL)
    memcpy(forward, a_begin, nseg * nt * sizeof(FP_TYPE));
  free(s_begin); free(s_end);
  free(a_begin); free(a_end);
  free(s_reseg); free(t_reseg);
  return reseg;
# undef sb
# undef se
# else
  if(ab_ != NULL)
    *ab_ = a_begin;
  else
    free(a_begin);
  return a_end;
# endif
}
