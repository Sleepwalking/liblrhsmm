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
  int* reseg   = calloc(nseg, sizeof(int));
#   define sb(t, i) s_begin[(t) * nseg + (i)]
#   define se(t, i) s_end  [(t) * nseg + (i)]
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
    pinit = p(t, j) + (j > 0 && t > 0 ? ae(t - 1, j - 1) : 0);
    if(j == 0 && t > 0) pinit = NEGINF;
    if(t == 0 && j > 0) pinit = NEGINF;
#   ifdef VITERBI
    maxfrom = j - 1;
#   endif
    
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
  reseg[nseg - 1] = se(nt - 1, nseg - 1);

  // backtracking
  int j = nseg - 1;
  while(j > 0) {
    int trans = sb(reseg[j], j);
    reseg[j - 1] = se(reseg[j] - 1, trans);
    j --;
  }
  
  // shift
  for(int i = 0; i < nseg - 1; i ++)
    reseg[i] = reseg[i + 1];
  reseg[nseg - 1] = nt;
# endif

  free(a_end);
# undef ab
# undef ae
# ifdef VITERBI
  if(forward != NULL)
    memcpy(forward, a_begin, nseg * nt * sizeof(FP_TYPE));
  free(s_begin); free(s_end);
  free(a_begin);
  return reseg;
#   undef sb
#   undef se
# else
  return a_begin;
# endif
}
