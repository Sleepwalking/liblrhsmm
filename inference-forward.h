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
  #define s(t, i) s_[(t) * nseg + (i)]
# endif

  for(int t = 0; t < nt; t ++)
  for(int i = 0; i < nseg; i ++) {
    a(t, i) = NEGINF;
# ifdef VITERBI
    s(t, i) = 0;
# endif
  }
  // t = 0
  a(0, 0) = p(0, 0);
  // t = 1, ..., nt-1
  for_tj_forward(1, 1, 1, 1)
    int jdur = seg -> time[j - 1] - (j == 1 ? 0 : seg -> time[j - 2]);
    int maxdur, mindur = 1;
    lrh_duration* srcdr = model -> durations[seg -> durstate[j - 1]];
    calculate_maxdur(srcdr, jdur);
    FP_TYPE lgsum = NEGINF;
#   ifdef VITERBI
    int maxfrom = max(0, t - jdur);
#   endif
    FP_TYPE pobserv = 0;

    // s_t = j, T_t, s_{t-d} = j-1, T_{t-d}
    int D = min(maxdur, t + 1);
    int d = mindur;
    if(d == 1 && d <= t) { // d = 1 case
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, 1);
      FP_TYPE plast = a(t - 1, j - 1);
#     ifdef VITERBI
      FP_TYPE pmul = pdur + plast; // pobserv = 0
      if(pmul > lgsum) {
        lgsum = pmul;
        maxfrom = t - 1;
      }
#     else
      lgsum = pdur + plast; // lgsum was zero before this line; pobserv = 0
#     endif
      d = 2;
    }
    for(; d < D; d ++) { // 1 < d <= t
      pobserv += p(t - d + 1, j - 1);
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, d);
      FP_TYPE plast = a(t - d, j - 1);
#     ifdef VITERBI
      FP_TYPE pmul = pobserv + pdur + plast;
      if(pmul > lgsum) {
        lgsum = pmul;
        maxfrom = t - d;
      }
#     else
      lgsum = lrh_lse(lgsum, pobserv + pdur + plast);
#     endif
    }
/*  unoptimized code (for reference)
    for(int d = mindur; d < maxdur; d ++) {
      if(d > t) break;
      if(d > 1)
        pobserv += p(t - d + 1, j - 1);
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, d);
      FP_TYPE plast = a(t - d, j - 1);
#     ifdef VITERBI
      FP_TYPE pmul = pobserv + pdur + plast;
      if(pmul > lgsum) {
        lgsum = pmul;
        maxfrom = t - d;
      }
#     else
      lgsum = lrh_lse(lgsum, pobserv + pdur + plast);
#     endif
    } */

    a(t, j) = lgsum + p(t, j);
#   ifdef VITERBI
    s(t, j) = maxfrom;
#   endif
  end_for_tj()

/*
//  if(lrh_debug_flag)
  for(int t = 0; t < nt; t ++) {
    FP_TYPE maxt = -9999999;
    for(int j = 0; j < nseg; j ++)
      maxt = max(p(t, j), maxt);
    printf("%d %f\n", t, maxt);
  }
*/

# ifdef VITERBI
  // Finalization
  int j = nseg - 1;
  int maxfrom = seg -> time[nseg - 2];
  FP_TYPE maxp = NEGINF;

  // maximize the total probability with regard to the onset of the last state
  FP_TYPE pobserv = 0;
  for(int d = 1; d < nt - 1; d ++) {
    pobserv += p(nt - d, nseg - 1);
    FP_TYPE p = pobserv + a(nt - d - 1, nseg - 1);
    if(p > maxp) {
      maxp = p;
      maxfrom = nt - d;
    }
  }
  reseg[nseg - 1] = maxfrom;

  // Backtracking
  while(j > 0) {
    reseg[j - 1] = s(reseg[j], j);
    j --;
  }
  
  // Shift
  for(int i = 0; i < nseg - 1; i ++)
    reseg[i] = reseg[i + 1];
  reseg[nseg - 1] = nt;
# endif

# ifdef VITERBI
  if(forward == NULL) free(a_);
  free(s_);
  return reseg;
# else
  return a_;
# endif
}

