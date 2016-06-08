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

#ifndef LRHSMM_DATA_H
#define LRHSMM_DATA_H

typedef struct {
  FP_TYPE** data;   // stream-wise array of matrices (time * dimension)
  int nstream, nt;
  int* ndim;       // dimension of each stream
} lrh_observ;

// access the content of a lrh_observ at time t, l-th stream, i-th dimension
#define lrh_obm(obj, t, i, l) (((obj) -> data[l])[(t) * ((obj) -> ndim[l]) + i])

typedef struct {
  int* time;       // right boundary of each segmentation
  int** outstate;  // stream-wise array of vectors of indices of output pdfs in lrh_stream
  int* durstate;   // vector of indices of duration pdfs in lrh_model
  int nstream;
  int nseg;
} lrh_seg;

// container for unsegmented data
typedef struct {
  lrh_observ** samples;
  int nsample;
} lrh_observset;

// container for segmentations
typedef struct {
  lrh_seg** samples;
  int nsample;
} lrh_segset;

typedef struct {
  lrh_observset* observset;
  lrh_segset* segset;
} lrh_dataset;

lrh_observ* lrh_create_observ(int nstream, int nt, int* ndim);
lrh_seg* lrh_create_seg(int nstream, int nseg);
lrh_seg* lrh_seg_copy(lrh_seg* src);
lrh_segset* lrh_segset_copy(lrh_segset* src);
void lrh_delete_observ(lrh_observ* dst);
void lrh_delete_seg(lrh_seg* dst);

lrh_observset* lrh_create_empty_observset(int nsample);
lrh_segset* lrh_create_empty_segset(int nsample);
void lrh_delete_observset(lrh_observset* dst);
void lrh_delete_segset(lrh_segset* dst);

// regroup observation and segmentation into smaller units
lrh_dataset* lrh_regroup_data(lrh_observset* srcobset, lrh_segset* srcsgset, int unitsize);

#endif

