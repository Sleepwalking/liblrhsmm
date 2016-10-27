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

#define p(t, i) p_[(t) * nseg + (i)]
#define a(t, i) a_[(t) * nseg + (i)]
#define b(t, i) b_[(t) * nseg + (i)]
#define x(t, i) x_[(t) * nseg + (i)]
#define y(t, i) y_[(t) * nseg + (i)]
#define m(t, i) m_[(t) * nseg + (i)]

#define alloc_init(d) ((float*)lrh_mempool_calloc(pool, d, sizeof(FP_TYPE)))
#define alloc(d) ((float*)lrh_mempool_malloc(pool, (d) * sizeof(FP_TYPE)))
#define x_safe(t, i, d) ((d > x(t, i).upper || d < x(t, i).lower) ? 0 : x(t, i).p[d - x(t, i).lower])

#define calculate_maxdur(srcdr, stdur) \
  maxdur = stdur * (1.0 + lrh_inference_duration_extra_factor) + lrh_inference_duration_extra; \
  maxdur = max(maxdur, srcdr -> mean + sqrt(srcdr -> var) * 2); \
  if(srcdr -> ceil > 0) maxdur = min(maxdur, srcdr -> ceil);

#define for_tj_forward(min_t, min_j, max_t, max_j) \
  for(int i = 0; i < nseg; i ++) { \
    int tbegin = i == 0 ? 0 : seg -> time[i - 1]; \
    for(int t = max(min_t, tbegin); t <= min(nt - max_t, seg -> time[i]); t ++) { \
      for(int j = max(min_j, i - lrh_inference_stprune); \
        j <= min(nseg - max_j, i + lrh_inference_stprune); j ++) {

#define for_tj_backward(min_t, min_j, max_t, max_j) \
  for(int i = nseg - 1; i >= 0; i --) { \
    int tbegin = i == 0 ? 0 : seg -> time[i - 1]; \
    for(int t = min(nt - max_t, seg -> time[i] - 1); t >= max(min_t, tbegin); t --) { \
      for(int j = max(min_j, i - lrh_inference_stprune); \
        j <= min(nseg - max_j, i + lrh_inference_stprune); j ++) {

#define end_for_tj() }}}

