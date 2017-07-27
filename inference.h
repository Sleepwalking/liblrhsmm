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

#ifndef LRHSMM_INFERENCE_H
#define LRHSMM_INFERENCE_H

#include "data.h"
#include "model.h"
#include "mempool.h"

typedef struct {
  float* p;
  int lower, upper;
} lrh_pslice;

FP_TYPE lrh_gmm_mixprob_lg(lrh_gmm* src, FP_TYPE* observ, int k);
double lrh_gmm_mixprob(lrh_gmm* src, FP_TYPE* observ, int k);
FP_TYPE lrh_gmm_outputprob_lg(lrh_gmm* src, FP_TYPE* observ);
double lrh_gmm_outputprob(lrh_gmm* src, FP_TYPE* observ);
FP_TYPE lrh_duration_prob_lg(lrh_duration* src, int duration);
double lrh_duration_prob(lrh_duration* src, FP_TYPE duration);

FP_TYPE* lrh_sample_outputprob_lg(lrh_model* model, lrh_observ* observ, lrh_seg* seg);
FP_TYPE* lrh_sample_outputprob_lg_full(lrh_model* model, lrh_observ* observ, lrh_seg* seg);
FP_TYPE* lrh_sample_outputprob(lrh_model* model, lrh_observ* observ, lrh_seg* seg);

/*
  ae: forward termination probability     p(s'_t = i, o_1, ..., o_t)
  ab: forward initiation probability      p(s*_t = i, o_1, ..., o_t)
  be: backward termination probability    p(o_{t+1}, ..., o_T | s'_t = i)
  bb: backward initiation probability     p(o_{t+1}, ..., o_T | s*_t = i)
  lrh_forward returns ae by default; ab is optional (to be allocated, which can be NULL)
  lrh_backward returns be by default; bb is optional (to be allocated, which can be NULL)
  In most cases ae/be instead of ab/bb are used for the following computation.
*/
FP_TYPE* lrh_forward(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt, FP_TYPE** ab);
FP_TYPE* lrh_forward_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt);
lrh_seg* lrh_viterbi(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt, FP_TYPE** ae);
int*     lrh_viterbi_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt, FP_TYPE** a);
FP_TYPE* lrh_backward(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt, FP_TYPE** bb);
FP_TYPE* lrh_backward_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt);

/*
  lrh_total: compute total probability from both forward and backward probabilities (HSMM)
  lrh_total_backward_b: compute total probability from backward termination probability (HSMM)
  lrh_total_pseudo: compute total probability from foward probability (HMM)
  lrh_total_pseudo_backward: compute total probability from backward probability (HMM)

  TL;DR: always use either lrh_total_backward_b or lrh_total evaluated at the first state to
    compute total probability.
  ---
  lrh_total_foward/lrh_total_forward_b/lrh_total_backward are not implemented because the
    definition of observation probability after the end of the sequence is not clear.
  For HSMMs, lrh_total_backward_b always gives accurate results.
  lrh_total gives accurate results when evaluated at the first state. However the assumption
    behind the total probability equation is violated when there are skipped states, giving
    wrong results after the first state. This is what the equation does. It's not a bug.
*/
FP_TYPE  lrh_total(FP_TYPE* a_, FP_TYPE* b_, int nt, int nseg, int i);
FP_TYPE  lrh_total_backward_b(FP_TYPE* b_, FP_TYPE* output_lg, int nseg);
FP_TYPE  lrh_total_pseudo(FP_TYPE* a_, int nt, int nseg);
FP_TYPE  lrh_total_pseudo_backward(FP_TYPE* b_, FP_TYPE* output_lg, int nseg);

/*
  lrh_durocp: compute consecutive state occupancy probability from ae and be probabilities (HSMM)
  lrh_stocp: compute single state occupancy probability from durocp (HSMM)
  lrh_mixocp: compute mixture occupancy probability from stocp (HSMM)
*/
lrh_pslice* lrh_durocp(lrh_model* model, lrh_seg* seg, FP_TYPE* ae, FP_TYPE* be,
  FP_TYPE* output_lg, FP_TYPE total, int nt, lrh_mempool* pool);
lrh_pslice* lrh_sttran_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* a, FP_TYPE* b,
  FP_TYPE* output_lg, FP_TYPE total, int nt, lrh_mempool* pool);
FP_TYPE   * lrh_stocp(lrh_model* model, lrh_seg* seg, lrh_pslice* durocp, int nt);
lrh_pslice* lrh_mixocp(lrh_model* model, lrh_observ* observ, lrh_seg* seg,
  FP_TYPE* stocp, int stream, lrh_mempool* pool);
FP_TYPE*    lrh_stocp_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* a, FP_TYPE* b,
  FP_TYPE total, int nt);
lrh_pslice* lrh_mixocp_geometric(lrh_model* model, lrh_observ* observ, lrh_seg* seg,
  FP_TYPE* stocp, int stream, lrh_mempool* pool);

#endif
