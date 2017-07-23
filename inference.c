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

#include "inference.h"
#include <string.h>
#include <stdlib.h>
#include "common.h"
#include "math-funcs.h"

#include "inference-helper.h"

FP_TYPE lrh_gmm_mixprob_lg(lrh_gmm* src, FP_TYPE* observ, int k) {
  FP_TYPE ex = 0;
  //FP_TYPE mul = 0;
  FP_TYPE *mean, *var;
  mean = & lrh_gmmu(src, k, 0);
  var = & lrh_gmmv(src, k, 0);
  for(int i = 0; i < src -> ndim; i ++) {
    ex += (observ[i] - mean[i]) * (observ[i] - mean[i]) / var[i];
    //mul -= fastlog(var[i] * 2.0 * M_PI);
  }
  FP_TYPE ret = src -> _tmp_term[k] * 0.5 - ex * 0.5;
  return isnan(ret) ? -999 : ret;
}

double lrh_gmm_mixprob(lrh_gmm* src, FP_TYPE* observ, int k) {
  return fastexp(lrh_gmm_mixprob_lg(src, observ, k));
}

FP_TYPE lrh_gmm_outputprob_lg(lrh_gmm* src, FP_TYPE* observ) {
  if(src -> nmix == 1) { // single mixture case
    return lrh_gmm_mixprob_lg(src, observ, 0);
  } else { // multi-mixture case
    double p = 0;
    for(int k = 0; k < src -> nmix; k ++) { // k: mixture number
      p += lrh_gmm_mixprob(src, observ, k) * src -> weight[k];
    }
    return max(fastlog(p), NEGINF);
  }
}

double lrh_gmm_outputprob(lrh_gmm* src, FP_TYPE* observ) {
  return fastexp(lrh_gmm_outputprob_lg(src, observ));
}

FP_TYPE lrh_duration_prob_lg(lrh_duration* src, int duration) {
  if(duration < lrh_precompute_duration)
    return src -> _tmp_prep[duration];
  return lrh_duration_logp(src, duration);
}

double lrh_duration_prob(lrh_duration* src, FP_TYPE duration) {
  return fastexp(lrh_duration_prob(src, duration));
}

FP_TYPE* lrh_sample_outputprob_lg(lrh_model* model, lrh_observ* observ, lrh_seg* seg) {
  lrh_model_precompute(model);

  int nt = observ -> nt;
  int nseg = seg -> nseg;
  FP_TYPE* p_ = calloc(sizeof(FP_TYPE), nt * nseg);
  memset(p_, 0, sizeof(FP_TYPE) * nt * nseg);
  for(int t = 0; t < nt; t ++)
    for(int i = 0; i < nseg; i ++)
      p(t, i) = NEGINF;

  for_tj_forward(0, 0, 1, 1)
    FP_TYPE pobserv = 0;
    for(int l = 0; l < observ -> nstream; l ++) {
      int jst = seg -> outstate[l][j];
      lrh_gmm* gmm = model -> streams[l] -> gmms[jst];
      pobserv += lrh_gmm_outputprob_lg(gmm, & lrh_obm(observ, t, 0, l)) *
        model -> streams[l] -> weight;
    }
    p(t, j) = pobserv;
  end_for_tj()

  return p_;
}

FP_TYPE* lrh_sample_outputprob_lg_full(lrh_model* model, lrh_observ* observ, lrh_seg* seg) {
  lrh_model_precompute(model);

  int nt = observ -> nt;
  int nseg = seg -> nseg;
  FP_TYPE* p_ = calloc(nt * nseg, sizeof(FP_TYPE));
  memset(p_, 0, sizeof(FP_TYPE) * nt * nseg);
  for(int t = 0; t < nt; t ++)
    for(int i = 0; i < nseg; i ++)
      p(t, i) = NEGINF;

  int prune_range = lrh_inference_stprune_full_slope * nseg;
  for(int t = 0; t < nt; t ++)
  for(int j = max(0, t * nseg / nt - prune_range);
          j < min(nseg, t * nseg / nt + prune_range); j ++) {
    FP_TYPE pobserv = 0;
    for(int l = 0; l < observ -> nstream; l ++) {
      int jst = seg -> outstate[l][j];
      lrh_gmm* gmm = model -> streams[l] -> gmms[jst];
      pobserv += lrh_gmm_outputprob_lg(gmm, & lrh_obm(observ, t, 0, l)) *
        model -> streams[l] -> weight;
    }
    p(t, j) = pobserv;
  }

  return p_;
}

FP_TYPE* lrh_sample_outputprob(lrh_model* model, lrh_observ* observ, lrh_seg* seg) {
  int nt = observ -> nt;
  int nseg = seg -> nseg;
  FP_TYPE* p_ = lrh_sample_outputprob_lg(model, observ, seg);
  for(int t = 0; t < nt; t ++)
    for(int i = 0; i < nseg; i ++)
      p(t, i) = fastexp(p(t, i));
  return p_;
}

FP_TYPE* lrh_forward(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt) {
  lrh_model_precompute(model);
  #undef VITERBI
  #include "inference-forward.h"
}

int* lrh_viterbi(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt,
  FP_TYPE** forward) {
  lrh_model_precompute(model);
  #define VITERBI
  #include "inference-forward.h"
}

FP_TYPE* lrh_forward_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt) {
  lrh_model_precompute(model);
  #undef VITERBI
  #include "inference-forward-geometric.h"
}

int* lrh_viterbi_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt,
  FP_TYPE** forward) {
  lrh_model_precompute(model);
  #define VITERBI
  #include "inference-forward-geometric.h"
}

FP_TYPE* lrh_backward(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt) {
  lrh_model_precompute(model);

  int nseg = seg -> nseg;
  FP_TYPE* p_ = output_lg;
  FP_TYPE* b_begin = calloc(nseg * nt, sizeof(FP_TYPE));
  FP_TYPE* b_end   = calloc(nseg * nt, sizeof(FP_TYPE));
# define bb(t, i) b_begin[(t) * nseg + (i)]
# define be(t, i) b_end  [(t) * nseg + (i)]

  for(int t = 0; t < nt; t ++)
    for(int i = 0; i < nseg; i ++) {
      bb(t, i) = NEGINF;
      be(t, i) = NEGINF;
    }

  for_tj_backward(0, 0, 1, 1)
    int jdur = seg -> time[j] - (j == 0 ? 0 : seg -> time[j - 1]);
    int maxdur, mindur = 1;
    lrh_duration* srcdr = model -> durations[seg -> durstate[j]];
    calculate_maxdur(srcdr, jdur);

    // compute termination probability
    FP_TYPE pterm = NEGINF;
    if(t + 1 < nt && j + 1 < nseg) {
      pterm = p(t + 1, j + 1) + bb(t + 1, j + 1);
    } else if(t + 1 == nt) {
      pterm = 0;
    }
    be(t, j) = pterm;

    // compute initiation probability
    FP_TYPE pinit = NEGINF;
    FP_TYPE pobserv = 0;
    if(j == nseg - 1) {
      pinit = 0;
      for(int d = 2; d <= nt - t; d ++)
        pinit += p(t + d - 1, j);
    } else {
      for(int d = 2; d < mindur; d ++) {
        int t1 = t + d - 1;
        pobserv += (t1 < nt) ? p(t1, j) : 0;
      }
      for(int d = mindur; d < maxdur; d ++) {
        int t1 = t + d - 1;
        pobserv += (t1 < nt && t1 > t) ? p(t1, j) : 0;
        FP_TYPE pdur = lrh_duration_prob_lg(srcdr, d);
        FP_TYPE porigin = t1 < nt ? be(t1, j) : 0;
        pinit = lrh_lse(pinit, pobserv + pdur + porigin);
      }
    }

    bb(t, j) = pinit;

  end_for_tj()

# undef bb
# undef be
  free(b_end);

  return b_begin;
}

FP_TYPE* lrh_backward_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt) {
  lrh_model_precompute(model);

  int nseg = seg -> nseg;
  FP_TYPE* p_ = output_lg;
  FP_TYPE* b_ = calloc(nt * nseg, sizeof(FP_TYPE));

  for(int t = 0; t < nt; t ++)
  for(int i = 0; i < nseg; i ++)
    b(t, i) = NEGINF;

  for(int i = 0; i < nseg; i ++)
    b(nt - 1, i) = 0;

  // t = nt-2, ..., 0
  int prune_range = lrh_inference_stprune_full_slope * nseg;
  for(int t = nt - 2; t >= 0; t --) {
    for(int j = max(0, t * nseg / nt - prune_range);
            j < min(nseg, t * nseg / nt + prune_range); j ++) {
      FP_TYPE pstay = 1.0 - 1.0 / max(1.01, model -> durations[seg -> durstate[j]] -> mean);
      
      FP_TYPE p = b(t + 1, j) + p(t + 1, j) + log(pstay);
      FP_TYPE lgsum = p;
      if(j != nseg - 1) {
        FP_TYPE ptrans = 1.0 / max(1.01, model -> durations[seg -> durstate[j]] -> mean);
        p = b(t + 1, j + 1) + p(t + 1, j + 1) + log(ptrans);
        lgsum = lrh_lse(lgsum, p);
      }

      b(t, j) = lgsum;
    }
  }
  return b_;
}

FP_TYPE lrh_total(FP_TYPE* a_, FP_TYPE* b_, int nt, int nseg, int i) {
  FP_TYPE lgsum = NEGINF;
  for(int t = 0; t < nt; t ++)
    lgsum = lrh_lse(lgsum, a(t, i) + b(t, i));
  return lgsum;
}

FP_TYPE lrh_total_forward(FP_TYPE* a_, FP_TYPE* output_lg, int nt, int nseg) {
  FP_TYPE lgsum = NEGINF;
  FP_TYPE* p_ = output_lg;
  FP_TYPE pobserv = 0;
  for(int t = nt - 1; t > nseg; t --) {
    pobserv += t == nt - 1 ? 0 : p(t + 1, nseg - 1);
    FP_TYPE pt = a(t, nseg - 1) + pobserv;
    lgsum = lrh_lse(lgsum, pt);
  }
  return lgsum;
}

FP_TYPE  lrh_total_backward(FP_TYPE* b_, FP_TYPE* output_lg, int nseg) {
  FP_TYPE* p_ = output_lg;
  return b(0, 0) + p(0, 0);
}

FP_TYPE lrh_total_pseudo(FP_TYPE* a_, int nt, int nseg) {
  FP_TYPE lgsum = NEGINF;
  for(int j = 0; j < nseg; j ++)
    lgsum = lrh_lse(lgsum, a(nt - 1, j));
  return lgsum;
}

FP_TYPE lrh_total_pseudo_backward(FP_TYPE* b_, FP_TYPE* output_lg, int nseg) {
  FP_TYPE* p_ = output_lg;
  return b(0, 0) + p(0, 0);
}

lrh_pslice* lrh_durocp(lrh_model* model, lrh_seg* seg, FP_TYPE* forward, FP_TYPE* backward,
  FP_TYPE* output_lg, FP_TYPE total, int nt, lrh_mempool* pool) {
  lrh_model_precompute(model);

  int nseg = seg -> nseg;
  FP_TYPE* p_ = output_lg;
  FP_TYPE* a_ = forward;
  FP_TYPE* b_ = backward;
  lrh_pslice* x_ = calloc(nt * nseg, sizeof(lrh_pslice));

  for(int t = 0; t < nt; t ++)
  for(int i = 0; i < nseg; i ++) {
    x(t, i).p = NULL;
    x(t, i).lower = 0;
    x(t, i).upper = -1;
  }

  // t = 0, ..., nt-1
  for_tj_forward(0, 0, 1, 2)
    int jdur = seg -> time[j] - (j == 0 ? 0 : seg -> time[j - 1]);
    int maxdur, mindur = 1;
    lrh_duration* srcdr = model -> durations[seg -> durstate[j]];
    calculate_maxdur(srcdr, jdur);
    if(j == nseg - 1) {
      mindur = nt - t;
      maxdur = mindur + lrh_inference_duration_extra;
    }

    x(t, j).p = alloc(maxdur - mindur);
    x(t, j).lower = mindur;
    x(t, j).upper = maxdur - 1;

    FP_TYPE pobserv = 0;
    int D = nt - t;
    int Dmid = min(maxdur, D);
    int d = mindur;
    if(d == 1) { // d = 1 case
      FP_TYPE pback = t + 1 < nt ? b(t + 1, j + 1) : 0;
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, 1);
      FP_TYPE pnext = t + 1 < nt ? p(t + 1, j + 1) : 0;
      FP_TYPE pd = a(t, j) + pback + pdur + pnext - total; // pobserv = 0
      x(t, j).p[0] = lrh_negexp(pd);
      d = 2;
    }
    for(; d < Dmid; d ++) { // d > 1 && t + d < nt
      pobserv += p(t + d - 1, j);
      FP_TYPE pback = b(t + d, j + 1);
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, d);
      FP_TYPE pnext = p(t + d, j + 1);
      FP_TYPE pd = a(t, j) + pback + pdur + pobserv + pnext - total;
      x(t, j).p[d - mindur] = lrh_negexp(pd);
    }
    if(d < maxdur && d == D) { // t + d = nt
      pobserv += p(t + d - 1, j);
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, d);
      FP_TYPE pd = a(t, j) + pdur + pobserv - total; // pback = pnext = 0
      x(t, j).p[d - mindur] = lrh_negexp(pd);
      d ++;
    }
    for(; d < maxdur; d ++) { // t + d > nt
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, d);
      FP_TYPE pd = a(t, j) + pdur + pobserv - total; // pback = pnext = 0
      x(t, j).p[d - mindur] = lrh_negexp(pd);
    }
    
    // trim x(t, j) to reduce inner loop in lrh_stocp
    int n = x(t, j).upper - x(t, j).lower;
    FP_TYPE p0 = x(t, j).p[0]; // to ensure the loop terminates at n = 0
    x(t, j).p[0] = 1; while(x(t, j).p[n] < 1e-5) n --; x(t, j).p[0] = p0;
    x(t, j).upper = n + x(t, j).lower;

/*  unoptimized code (for reference)
    for(int d = mindur; d < maxdur; d ++) {
      if(d > 1 && t + d <= nt)
        pobserv += p(t + d - 1, j);
      FP_TYPE pback = t + d < nt ? b(t + d, j + 1) : 0;
      FP_TYPE pdur = lrh_duration_prob_lg(srcdr, d);
      FP_TYPE pnext = t + d < nt ? p(t + d, j + 1) : 0;
      FP_TYPE pd = a(t, j) + pback + pdur + pobserv + pnext - total;
      x(t, j).p[d - mindur] = fastexp(min(0.0, pd));
    } */
  end_for_tj()

  return x_;
}

FP_TYPE* lrh_stocp(lrh_model* model, lrh_seg* seg, lrh_pslice* durocp, int nt) {
  lrh_model_precompute(model);

  int nseg = seg -> nseg;
  lrh_pslice* x_ = durocp;
  FP_TYPE* y_ = calloc(nt * nseg, sizeof(FP_TYPE));

  for(int t = 0; t < nt; t ++)
    for(int i = 0; i < nseg; i ++)
      y(t, i) = 0;

  // t = 0, ..., nt-1
  for_tj_forward(0, 0, 1, 1)
    int jdur = seg -> time[j] - (j == 0 ? 0 : seg -> time[j - 1]);
    int maxdur;
    lrh_duration* srcdr = model -> durations[seg -> durstate[j]];
    calculate_maxdur(srcdr, jdur);
    FP_TYPE xsum = 0;
    if(j == 0) {
      for(int d = 1 + t; d < maxdur + t; d ++)
        xsum += x_safe(0, 0, d);
    } else {
      for(int u = 0; u < min(maxdur, t); u ++) {
        lrh_pslice slice = x(t - u, j);
        int base = u - slice.lower;
        int vl = max(1, -base), vh = min(maxdur - 1, slice.upper - u);
        for(int v = vl; v <= vh; v ++)
          xsum += slice.p[v + base];
      }
    }
    y(t, j) = min(1.0, xsum);
  end_for_tj()

  return y_;
}

FP_TYPE* lrh_stocp_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* a_, FP_TYPE* b_,
  FP_TYPE total, int nt) {
  lrh_model_precompute(model);

  int nseg = seg -> nseg;
  FP_TYPE* y_ = calloc(nt * nseg, sizeof(FP_TYPE));

  for(int t = 0; t < nt; t ++)
    for(int i = 0; i < nseg; i ++)
      y(t, i) = 0;

  // t = 0, ..., nt-1
  int prune_range = lrh_inference_stprune_full_slope * nseg;
  for(int t = 1; t < nt - 1; t ++) {
    for(int j = max(0, t * nseg / nt - prune_range);
            j < min(nseg, t * nseg / nt + prune_range); j ++) {
      FP_TYPE ly = a(t, j) + b(t, j) - total;
      y(t, j) = lrh_negexp(ly);
      y(t, j) = y(t, j);
    }
  }

  return y_;
}

lrh_pslice* lrh_sttran_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* a_, FP_TYPE* b_,
  FP_TYPE* p_, FP_TYPE total, int nt, lrh_mempool* pool) {
  lrh_model_precompute(model);

  int nseg = seg -> nseg;
  lrh_pslice* y_ = calloc(nt * nseg, sizeof(lrh_pslice));

  // t = 0, ..., nt-1
  int prune_range = lrh_inference_stprune_full_slope * nseg;
  for(int t = 1; t < nt - 1; t ++) {
    for(int j = max(0, t * nseg / nt - prune_range);
            j < min(nseg, t * nseg / nt + prune_range); j ++) {
      FP_TYPE pstay = 1.0 - 1.0 / max(1.01, model -> durations[seg -> durstate[j]] -> mean);
      y(t, j).p = alloc(2);
      y(t, j).upper = 2;
      y(t, j).p[0] = lrh_negexp(a(t - 1, j) + b(t, j) + p(t, j) + log(pstay) - total);
      y(t, j).p[1] = 1;
      if(j > 0) {
        FP_TYPE ptrans = 1.0 / max(1.01, model -> durations[seg -> durstate[j - 1]] -> mean);
        y(t, j).p[1] = lrh_negexp(a(t - 1, j - 1) + b(t, j) + p(t, j) + log(ptrans) - total);
      }
    }
  }

  return y_;
}

lrh_pslice* lrh_mixocp(lrh_model* model, lrh_observ* observ, lrh_seg* seg,
  FP_TYPE* stocp, int stream, lrh_mempool* pool) {
  lrh_model_precompute(model);

  int nt = observ -> nt;
  int nseg = seg -> nseg;
  FP_TYPE* y_ = stocp;
  lrh_pslice* m_ = calloc(nt * nseg, sizeof(lrh_pslice));

  for(int t = 0; t < nt; t ++)
  for(int i = 0; i < nseg; i ++) {
    m(t, i).p = NULL;
    m(t, i).lower = 0;
    m(t, i).upper = 0;
  }

  for_tj_forward(0, 0, 1, 1)
    lrh_gmm* srcmx = model -> streams[stream] -> gmms[seg -> outstate[stream][j]];
    m(t, j).p = alloc(srcmx -> nmix);
    m(t, j).upper = srcmx -> nmix;
    FP_TYPE mop = fastlog(y(t, j)) - lrh_gmm_outputprob_lg(srcmx, & lrh_obm(observ, t, 0, stream));
    for(int k = 0; k < srcmx -> nmix; k ++) {
      m(t, j).p[k] = srcmx -> weight[k] *
        lrh_exp(mop + lrh_gmm_mixprob_lg(srcmx, & lrh_obm(observ, t, 0, stream), k));
      if(isnan(m(t, j).p[k])) // in case that the probability is too low
        m(t, j).p[k] = 0;
    }
  end_for_tj()

  return m_;
}

lrh_pslice* lrh_mixocp_geometric(lrh_model* model, lrh_observ* observ, lrh_seg* seg,
  FP_TYPE* stocp, int stream, lrh_mempool* pool) {
  lrh_model_precompute(model);

  int nt = observ -> nt;
  int nseg = seg -> nseg;
  FP_TYPE* y_ = stocp;
  lrh_pslice* m_ = calloc(nt * nseg, sizeof(lrh_pslice));

  for(int t = 0; t < nt; t ++)
  for(int i = 0; i < nseg; i ++) {
    m(t, i).p = NULL;
    m(t, i).lower = 0;
    m(t, i).upper = 0;
  }

  int prune_range = lrh_inference_stprune_full_slope * nseg;
  for(int t = 1; t < nt - 1; t ++) {
    for(int j = max(0, t * nseg / nt - prune_range);
            j < min(nseg, t * nseg / nt + prune_range); j ++) {
      lrh_gmm* srcmx = model -> streams[stream] -> gmms[seg -> outstate[stream][j]];
      m(t, j).p = alloc(srcmx -> nmix);
      m(t, j).upper = srcmx -> nmix;
      FP_TYPE mop = fastlog(y(t, j)) - lrh_gmm_outputprob_lg(srcmx, & lrh_obm(observ, t, 0, stream));
      for(int k = 0; k < srcmx -> nmix; k ++) {
        m(t, j).p[k] = srcmx -> weight[k] *
          lrh_negexp(mop + lrh_gmm_mixprob_lg(srcmx, & lrh_obm(observ, t, 0, stream), k));
        if(isnan(m(t, j).p[k])) // in case that the probability is too low
          m(t, j).p[k] = 0;
      }
    }
  }

  return m_;
}
