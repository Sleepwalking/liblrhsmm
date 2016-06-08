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

#include "common.h"

FP_TYPE lrh_duration_vfloor = 0.05;
FP_TYPE lrh_output_vfloor = 0.01;
FP_TYPE lrh_inference_duration_weight = 1.0;
FP_TYPE lrh_inference_stprune = 5;
FP_TYPE lrh_inference_stprune_full_slope = 0.3;
int     lrh_inference_duration_extra = 30;
FP_TYPE lrh_inference_duration_extra_factor = 1;
int     lrh_precompute_duration = 300;

