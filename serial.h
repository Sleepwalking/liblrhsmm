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

#ifndef LRHSMM_SERIAL_H
#define LRHSMM_SERIAL_H

// data serialization & deserialization

#include <stdint.h>
#include <stdlib.h>
typedef char bool;

#include "model.h"
#include "data.h"
#include "external/cmp/cmp.h"

int lrh_write_duration(cmp_ctx_t* dst, lrh_duration* src);
int lrh_write_gmm(cmp_ctx_t* dst, lrh_gmm* src);
int lrh_write_stream(cmp_ctx_t* dst, lrh_stream* src);
int lrh_write_model(cmp_ctx_t* dst, lrh_model* src);

lrh_duration* lrh_read_duration(cmp_ctx_t* src);
lrh_gmm* lrh_read_gmm(cmp_ctx_t* src);
lrh_stream* lrh_read_stream(cmp_ctx_t* src);
lrh_model* lrh_read_model(cmp_ctx_t* src);

int lrh_write_observ(cmp_ctx_t* dst, lrh_observ* src);
int lrh_write_observset(cmp_ctx_t* dst, lrh_observset* src);
int lrh_write_seg(cmp_ctx_t* dst, lrh_seg* src);
int lrh_write_segset(cmp_ctx_t* dst, lrh_segset* src);

lrh_observ* lrh_read_observ(cmp_ctx_t* src);
lrh_observset* lrh_read_observset(cmp_ctx_t* src);
lrh_seg* lrh_read_seg(cmp_ctx_t* src);
lrh_segset* lrh_read_segset(cmp_ctx_t* src);

#endif

