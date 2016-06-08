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

#ifndef LRHSMM_INFERENCE_H
#define LRHSMM_INFERENCE_H

#include "data.h"
#include "model.h"
#include "mempool.h"

typedef struct {
  float* p;
  int lower, upper;
} lrh_pslice;

FP_TYPE* lrh_sample_outputprob_lg(lrh_model* model, lrh_observ* observ, lrh_seg* seg);
FP_TYPE* lrh_sample_outputprob_lg_full(lrh_model* model, lrh_observ* observ, lrh_seg* seg);
FP_TYPE* lrh_sample_outputprob(lrh_model* model, lrh_observ* observ, lrh_seg* seg);

FP_TYPE* lrh_forward(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt);
FP_TYPE* lrh_forward_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt);
int*     lrh_viterbi(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt, FP_TYPE** forward);
int*     lrh_viterbi_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt, FP_TYPE** forward);
FP_TYPE* lrh_backward(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt);
FP_TYPE* lrh_backward_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* output_lg, int nt);

FP_TYPE  lrh_total(FP_TYPE* a_, FP_TYPE* b_, int nt, int nseg, int i);
FP_TYPE  lrh_total_forward(FP_TYPE* a_, FP_TYPE* output_lg, int nt, int nseg);
FP_TYPE  lrh_total_backward(FP_TYPE* b_, FP_TYPE* output_lg, int nseg);
FP_TYPE  lrh_total_pseudo(FP_TYPE* a_, int nt, int nseg);
FP_TYPE  lrh_total_pseudo_backward(FP_TYPE* b_, FP_TYPE* output_lg, int nseg);

lrh_pslice* lrh_durocp(lrh_model* model, lrh_seg* seg, FP_TYPE* forward, FP_TYPE* backward,
  FP_TYPE* output_lg, FP_TYPE total, int nt, lrh_mempool* pool);
lrh_pslice* lrh_sttran_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* forward, FP_TYPE* backward,
  FP_TYPE* output_lg, FP_TYPE total, int nt, lrh_mempool* pool);
FP_TYPE   * lrh_stocp(lrh_model* model, lrh_seg* seg, lrh_pslice* durocp, int nt);
lrh_pslice* lrh_mixocp(lrh_model* model, lrh_observ* observ, lrh_seg* seg,
  FP_TYPE* stocp, int stream, lrh_mempool* pool);
FP_TYPE*    lrh_stocp_geometric(lrh_model* model, lrh_seg* seg, FP_TYPE* forward, FP_TYPE* backward,
  FP_TYPE total, int nt);
lrh_pslice* lrh_mixocp_geometric(lrh_model* model, lrh_observ* observ, lrh_seg* seg,
  FP_TYPE* stocp, int stream, lrh_mempool* pool);

#endif

