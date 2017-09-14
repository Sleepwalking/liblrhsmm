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

#include <stdio.h>
#include "estimate.h"
#include "common.h"
#include "math-funcs.h"

#include "inference-helper.h"

void lrh_duration_update(lrh_duration* dst, lrh_duration_stat* src, int geometric) {
  if(src -> ocpsum < EPS) return;
  if(geometric) {
    dst -> mean = fmax(1.0, 1.0 / (src -> valsum / src -> ocpsum));
    dst -> var = dst -> mean;
  } else {
    dst -> mean = src -> valsum / src -> ocpsum;
    if(dst -> fixed_mean >= 0) dst -> mean = dst -> fixed_mean;
    dst -> var = src -> sqrsum / src -> ocpsum - dst -> mean * dst -> mean;
  }
  dst -> var = fmax(dst -> var, lrh_duration_vfloor);
  dst -> var = fmax(dst -> var, dst -> vfloor);
  dst -> _tmp_term = LRH_UNINITIALIZED;
}

void lrh_gmm_update(lrh_gmm* dst, lrh_gmm_stat* src, int geometric) {
  if(src -> ocpsum < EPS) return;
  for(int k = 0; k < src -> nmix; k ++) {
    FP_TYPE wk = src -> weightsum[k];
    if(wk < EPS) continue;
    dst -> weight[k] = wk / src -> ocpsum;
    for(int i = 0; i < src -> ndim; i ++) {
      lrh_gmmu(dst, k, i) = lrh_gmmvs(src, k, i) / wk;
      lrh_gmmv(dst, k, i) = lrh_gmmss(src, k, i) / wk - lrh_gmmu(dst, k, i) * lrh_gmmu(dst, k, i);
      lrh_gmmv(dst, k, i) = fmax(lrh_gmmv(dst, k, i), lrh_output_vfloor);
      lrh_gmmv(dst, k, i) = fmax(lrh_gmmv(dst, k, i), lrh_gmmvf(dst, k, i));
    }
    dst -> _tmp_term[k] = LRH_UNINITIALIZED;
  }
}

void lrh_stream_update(lrh_stream* dst, lrh_stream_stat* src, int geometric) {
  for(int i = 0; i < dst -> ngmm; i ++)
    lrh_gmm_update(dst -> gmms[i], src -> gmms[i], geometric);
}

void lrh_model_update(lrh_model* dst, lrh_model_stat* src, int geometric) {
  for(int i = 0; i < dst -> nstream; i ++)
    lrh_stream_update(dst -> streams[i], src -> streams[i], geometric);
  for(int i = 0; i < dst -> nduration; i ++)
    lrh_duration_update(dst -> durations[i], src -> durations[i], geometric);
}

void lrh_collect_init(lrh_model_stat* dst, lrh_observ* srcob, lrh_seg* srcseg) {
  for(int s = 0; s < srcseg -> nseg; s ++) {
    int tbegin = s == 0 ? 0 : srcseg -> time[s - 1];
    int tend = srcseg -> time[s];
    lrh_duration_stat* dstdr = dst -> durations[srcseg -> durstate[s]];
    dstdr -> valsum += tend - tbegin;
    dstdr -> sqrsum += (tend - tbegin) * (tend - tbegin);
    dstdr -> ocpsum ++;
    for(int t = tbegin; t < tend; t ++)
      for(int l = 0; l < dst -> nstream; l ++) {
        lrh_gmm_stat* dstgmm = dst -> streams[l] -> gmms[srcseg -> outstate[l][s]];
        if(dstgmm -> nmix != 1) {
          fprintf(stderr, "Error: lrh_model_init doesn't support multi-mixture GMM "
            "at the current stage.\n");
          return;
        }
        dstgmm -> ocpsum ++;
        dstgmm -> weightsum[0] ++;
        for(int i = 0; i < srcob -> ndim[l]; i ++) {
          lrh_gmmvs(dstgmm, 0, i) += lrh_obm(srcob, t, i, l);
          lrh_gmmss(dstgmm, 0, i) += lrh_obm(srcob, t, i, l) * lrh_obm(srcob, t, i, l);
        }
      }
  }
}

void lrh_collect_ocp(lrh_model_stat* dst, lrh_observ* observ, lrh_seg* seg,
  lrh_pslice* durocp, FP_TYPE* stocp, lrh_pslice** mixocp) {
  int nt = observ -> nt;
  int nseg = seg -> nseg;
  lrh_pslice* x_ = durocp;
  FP_TYPE*     y_ = stocp;

  for_tj_forward(0, 0, 1, 1)
    lrh_duration_stat* dstdr = dst -> durations[seg -> durstate[j]];
    // duration
    for(int d = x(t, j).lower; d <= x(t, j).upper; d ++) {
      FP_TYPE pd = x(t, j).p[d - x(t, j).lower];
      dstdr -> valsum += pd * d;
      dstdr -> sqrsum += pd * d * d;
      dstdr -> ocpsum += pd;
    }

    for(int l = 0; l < dst -> nstream; l ++) {
      lrh_pslice* m_ = mixocp[l];
      int jst = seg -> outstate[l][j];
      lrh_gmm_stat* dstmx = dst -> streams[l] -> gmms[jst];

      dstmx -> ocpsum += y(t, j);

      // mixtures
      if(m(t, j).p != NULL)
      for(int k = 0; k < dstmx -> nmix; k ++) {
        FP_TYPE pk = m(t, j).p[k];
        FP_TYPE* o = & lrh_obm(observ, t, 0, l);
        dstmx -> weightsum[k] += pk;
        for(int n = 0; n < dstmx -> ndim; n ++) {
          lrh_gmmvs(dstmx, k, n) += pk * o[n];
          lrh_gmmss(dstmx, k, n) += pk * o[n] * o[n];
        }
      }
    }

  end_for_tj()
}

void lrh_collect_ocp_geometric(lrh_model_stat* dst, lrh_observ* observ, lrh_seg* seg,
  lrh_pslice* sttran, FP_TYPE* stocp, lrh_pslice** mixocp) {
  int nt = observ -> nt;
  int nseg = seg -> nseg;
  lrh_pslice* x_ = sttran;
  FP_TYPE* y_ = stocp;

  int prune_range = lrh_inference_stprune_full_slope * nseg;
  for(int t = 1; t < nt - 1; t ++) {
    for(int j = max(0, t * nseg / nt - prune_range);
            j < min(nseg, t * nseg / nt + prune_range); j ++) {
      lrh_duration_stat* dstdr = dst -> durations[seg -> durstate[j]];
      if(x(t, j).p != NULL)
        dstdr -> valsum += x(t, j).p[1]; // used to store ptrans
      dstdr -> ocpsum += y(t, j);

      for(int l = 0; l < dst -> nstream; l ++) {
        lrh_pslice* m_ = mixocp[l];
        int jst = seg -> outstate[l][j];
        lrh_gmm_stat* dstmx = dst -> streams[l] -> gmms[jst];

        dstmx -> ocpsum += y(t, j);

        // mixtures
        if(m(t, j).p != NULL)
        for(int k = 0; k < dstmx -> nmix; k ++) {
          FP_TYPE pk = m(t, j).p[k];
          FP_TYPE* o = & lrh_obm(observ, t, 0, l);
          dstmx -> weightsum[k] += pk;
          for(int n = 0; n < dstmx -> ndim; n ++) {
            lrh_gmmvs(dstmx, k, n) += pk * o[n];
            lrh_gmmss(dstmx, k, n) += pk * o[n] * o[n];
          }
        }
      }
    }
  }
}

FP_TYPE lrh_estimate(lrh_model_stat* dst, lrh_model* h, lrh_observ* observ, lrh_seg* seg) {
  lrh_mempool* pool = lrh_create_mempool(1024 * 1024);
  FP_TYPE* outp     = lrh_sample_outputprob_lg(h, observ, seg);
  FP_TYPE* a        = lrh_forward(h, seg, outp, observ -> nt, NULL);
  FP_TYPE* b        = lrh_backward(h, seg, outp, observ -> nt, NULL);
  FP_TYPE total     = lrh_total(a, b, observ -> nt, seg -> nseg, 0);

  lrh_pslice* docp  = lrh_durocp(h, seg, a, b, outp, total, observ -> nt, pool);
  free(a); free(b); free(outp);
  FP_TYPE* socp     = lrh_stocp(h, seg, docp, observ -> nt);

  lrh_pslice** mocp = calloc(h -> nstream, sizeof(lrh_pslice*));
  for(int i = 0; i < h -> nstream; i ++)
    mocp[i] = lrh_mixocp(h, observ, seg, socp, i, pool);

# pragma omp critical
  lrh_collect_ocp(dst, observ, seg, docp, socp, mocp);
  
  for(int i = 0; i < h -> nstream; i ++)
    free(mocp[i]);
  free(docp); free(socp); free(mocp);
  lrh_delete_mempool(pool);
  return total;
}

FP_TYPE lrh_estimate_geometric(lrh_model_stat* dst, lrh_model* h, lrh_observ* observ, lrh_seg* seg) {
  lrh_mempool* pool = lrh_create_mempool(1024 * 1024);
  FP_TYPE* outp     = lrh_sample_outputprob_lg_full(h, observ, seg);
  FP_TYPE* a        = lrh_forward_geometric(h, seg, outp, observ -> nt);
  FP_TYPE* b        = lrh_backward_geometric(h, seg, outp, observ -> nt);
  FP_TYPE total     = lrh_total_pseudo_backward(b, outp, seg -> nseg);
  lrh_pslice* strp  = lrh_sttran_geometric(h, seg, a, b, outp, total, observ -> nt, pool);
  free(outp);
  FP_TYPE* socp     = lrh_stocp_geometric(h, seg, a, b, total, observ -> nt);
  free(a); free(b);

  lrh_pslice** mocp = calloc(h -> nstream, sizeof(lrh_pslice*));
  for(int i = 0; i < h -> nstream; i ++)
    mocp[i] = lrh_mixocp_geometric(h, observ, seg, socp, i, pool);

# pragma omp critical
  lrh_collect_ocp_geometric(dst, observ, seg, strp, socp, mocp);
  
  for(int i = 0; i < h -> nstream; i ++)
    free(mocp[i]);
  free(strp); free(socp); free(mocp);
  lrh_delete_mempool(pool);
  return total;
}

