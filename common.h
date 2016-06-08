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

#ifndef LIBLRHSMM_COMMON_H
#define LIBLRHSMM_COMMON_H

#define LRH_UNINITIALIZED ((FP_TYPE)(-1e16))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

extern FP_TYPE lrh_duration_vfloor;
extern FP_TYPE lrh_output_vfloor;
extern FP_TYPE lrh_inference_duration_weight;
extern FP_TYPE lrh_inference_stprune;
extern FP_TYPE lrh_inference_stprune_full_slope;
extern int     lrh_inference_duration_extra;
extern FP_TYPE lrh_inference_duration_extra_factor;
extern int     lrh_precompute_duration;

#endif

