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

#ifndef LRHSMM_ESTIMATE_H
#define LRHSMM_ESTIMATE_H

#include "data.h"
#include "model.h"
#include "inference.h"

// functions for updating model parameters in maximization step
void lrh_duration_update(lrh_duration* dst, lrh_duration_stat* src, int geometric);
void lrh_gmm_update(lrh_gmm* dst, lrh_gmm_stat* src, int geometric);
void lrh_stream_update(lrh_stream* dst, lrh_stream_stat* src, int geometric);
void lrh_model_update(lrh_model* dst, lrh_model_stat* src, int geometric);

// collect expectation from a pair of observation and segmentation
void lrh_collect_init(lrh_model_stat* dst, lrh_observ* srcob, lrh_seg* srcseg);

// collect expectation from inference results
void lrh_collect_ocp(lrh_model_stat* dst, lrh_observ* observ, lrh_seg* seg,
  lrh_pslice* durocp, FP_TYPE* stocp, lrh_pslice** mixocp);
void lrh_collect_ocp_geometric(lrh_model_stat* dst, lrh_observ* observ, lrh_seg* seg,
  lrh_pslice* sttran, FP_TYPE* stocp, lrh_pslice** mixocp);

// run all relevant inference algorithms and then collect expectation from inference results
FP_TYPE lrh_estimate(lrh_model_stat* dst, lrh_model* h, lrh_observ* observ, lrh_seg* seg);
FP_TYPE lrh_estimate_geometric(lrh_model_stat* dst, lrh_model* h, lrh_observ* observ, lrh_seg* seg);

#endif

