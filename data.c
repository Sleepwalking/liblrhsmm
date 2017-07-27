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

#include <string.h>
#include <stdlib.h>
#include "data.h"

static void resize_arr(int** int_arr, FP_TYPE** fp_arr, int size) {
  *int_arr = realloc(*int_arr, size * sizeof(int));
  *fp_arr  = realloc(*fp_arr , size * sizeof(FP_TYPE));
}

static int getsize_arr(int* int_arr, int termination) {
  int size = 0;
  do {
    size ++;
    int_arr ++;
  } while(int_arr[-1] != termination);
  return size;
}

static void copy_arr(int* dst_int, FP_TYPE* dst_fp, int* src_int, FP_TYPE* src_fp,
  int termination) {
  do {
    *dst_int = *src_int;
    *dst_fp  = *src_fp;
    dst_int ++;
    dst_fp  ++;
    src_int ++;
    src_fp  ++;
  } while(src_int[-1] != termination);
}

static void append_arr(int** dst_int, FP_TYPE** dst_fp, int new_int, FP_TYPE new_fp,
  int termination) {
  int oldsize = getsize_arr(*dst_int, termination);
  if(new_int == termination) {
    dst_fp[0][oldsize - 1] = new_fp;
    return;
  }
  resize_arr(dst_int, dst_fp, oldsize + 1);
  dst_int[0][oldsize] = dst_int[0][oldsize - 1];
  dst_fp [0][oldsize] = dst_fp [0][oldsize - 1];
  dst_int[0][oldsize - 1] = new_int;
  dst_fp [0][oldsize - 1] = new_fp;
}

lrh_observ* lrh_create_observ(int nstream, int nt, int* ndim) {
  lrh_observ* ret = malloc(sizeof(lrh_observ));
  ret -> data = calloc(nstream, sizeof(FP_TYPE*));
  ret -> nstream = nstream;
  ret -> nt = nt;
  ret -> ndim = calloc(nstream, sizeof(int));
  memcpy(ret -> ndim, ndim, nstream * sizeof(int));
  for(int i = 0; i < nstream; i ++)
    ret -> data[i] = calloc(ndim[i] * nt, sizeof(FP_TYPE));
  return ret;
}

lrh_seg* lrh_create_seg(int nstream, int nseg) {
  lrh_seg* ret = malloc(sizeof(lrh_seg));
  ret -> time  = calloc(nseg, sizeof(int));
  ret -> outstate = calloc(nstream, sizeof(int*));
  ret -> durstate = calloc(nseg, sizeof(int));
  ret -> nseg = nseg;
  ret -> nstream = nstream;
  for(int i = 0; i < nstream; i ++)
    ret -> outstate[i] = calloc(nseg, sizeof(int));

  ret -> djump_out = calloc(nseg, sizeof(int*));
  ret -> pjump_out = calloc(nseg, sizeof(FP_TYPE*));
  ret -> djump_in  = calloc(nseg, sizeof(int*));
  ret -> pjump_in  = calloc(nseg, sizeof(FP_TYPE*));
  for(int i = 0; i < nseg; i ++) {
    ret -> djump_out[i] = malloc(sizeof(int));
    ret -> djump_in [i] = malloc(sizeof(int));
    ret -> pjump_out[i] = malloc(sizeof(FP_TYPE));
    ret -> pjump_in [i] = malloc(sizeof(FP_TYPE));
    ret -> djump_out[i][0] = 1;
    ret -> pjump_out[i][0] = 1.0;
    ret -> djump_in [i][0] = -1;
    ret -> pjump_in [i][0] = 1.0;
  }
  return ret;
}

void lrh_delete_observ(lrh_observ* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> nstream; i ++)
    if(dst -> data[i] != NULL)
      free(dst -> data[i]);
  free(dst -> data);
  free(dst -> ndim);
  free(dst);
}

void lrh_delete_seg(lrh_seg* dst) {
  if(dst == NULL) return;
  if(dst -> outstate != NULL)
    for(int i = 0; i < dst -> nstream; i ++)
      if(dst -> outstate[i] != NULL)
        free(dst -> outstate[i]);
  free(dst -> time);
  free(dst -> outstate);
  free(dst -> durstate);

  if(dst -> djump_out != NULL)
  for(int i = 0; i < dst -> nseg; i ++) {
    free(dst -> djump_out[i]);
    free(dst -> pjump_out[i]);
  }
  if(dst -> djump_in != NULL)
  for(int i = 0; i < dst -> nseg; i ++) {
    free(dst -> djump_in[i]);
    free(dst -> pjump_in[i]);
  }
  free(dst -> djump_out);
  free(dst -> pjump_out);
  free(dst);
}

lrh_observset* lrh_create_empty_observset(int nsample) {
  lrh_observset* ret = malloc(sizeof(lrh_observset));
  ret -> samples = calloc(nsample, sizeof(lrh_observ*));
  ret -> nsample = nsample;
  return ret;
}

lrh_segset* lrh_create_empty_segset(int nsample) {
  lrh_segset* ret = malloc(sizeof(lrh_segset));
  ret -> samples = calloc(nsample, sizeof(lrh_seg*));
  ret -> nsample = nsample;
  return ret;
}

lrh_seg* lrh_seg_copy(lrh_seg* src) {
  lrh_seg* ret = lrh_create_seg(src -> nstream, src -> nseg);
  for(int i = 0; i < src -> nseg; i ++) {
    ret -> time[i] = src -> time[i];
    ret -> durstate[i] = src -> durstate[i];
    for(int k = 0; k < src -> nstream; k ++)
      ret -> outstate[k][i] = src -> outstate[k][i];
    int n_out = getsize_arr(src -> djump_out[i], 1);
    int n_in  = getsize_arr(src -> djump_in [i], -1);
    resize_arr(& ret -> djump_out[i], & ret -> pjump_out[i], n_out);
    resize_arr(& ret -> djump_in [i], & ret -> pjump_in [i], n_in);
    copy_arr(ret -> djump_out[i], ret -> pjump_out[i],
             src -> djump_out[i], src -> pjump_out[i], 1);
    copy_arr(ret -> djump_in [i], ret -> pjump_in [i],
             src -> djump_in [i], src -> pjump_in [i], -1);
  }
  return ret;
}

lrh_segset* lrh_segset_copy(lrh_segset* src) {
  lrh_segset* ret = lrh_create_empty_segset(src -> nsample);
  for(int i = 0; i < src -> nsample; i ++)
    ret -> samples[i] = lrh_seg_copy(src -> samples[i]);
  return ret;
}

void lrh_seg_buildjumps(lrh_seg* dst) {
  for(int i = 0; i < dst -> nseg; i ++) { // reset in-bound transitions
    dst -> djump_in[i][0] = -1;
    dst -> pjump_in[i][0] = 1.0;
  }
  for(int i = 0; i < dst -> nseg; i ++) {
    int j = 0;
    do {
      int d = dst -> djump_out[i][j];
      int k = i + d;
      FP_TYPE p = dst -> pjump_out[i][j];
      // add i -> k tranition to k
      if(k >= 0 && k < dst -> nseg)
        append_arr(& dst -> djump_in[k], & dst -> pjump_in[k], -d, p, -1);
      j ++;
    } while(dst -> djump_out[i][j - 1] != 1);
  }
}

void lrh_delete_observset(lrh_observset* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> nsample; i ++)
    lrh_delete_observ(dst -> samples[i]);
  free(dst -> samples);
  free(dst);
}

void lrh_delete_segset(lrh_segset* dst)  {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> nsample; i ++)
    lrh_delete_seg(dst -> samples[i]);
  free(dst -> samples);
  free(dst);
}

lrh_dataset* lrh_regroup_data(lrh_observset* srcobset, lrh_segset* srcsgset, int unitsize) {
  if(srcobset -> nsample != srcsgset -> nsample) return NULL;
  int nsample = 0;
  for(int i = 0; i < srcsgset -> nsample; i ++)
    nsample += srcsgset -> samples[i] -> nseg / unitsize;
  lrh_dataset* ret = malloc(sizeof(lrh_dataset));
  ret -> observset = lrh_create_empty_observset(nsample);
  ret -> segset = lrh_create_empty_segset(nsample);

  nsample = 0;
  for(int i = 0; i < srcsgset -> nsample; i ++) {
    lrh_seg* srcsg = srcsgset -> samples[i];
    lrh_observ* srcob = srcobset -> samples[i];
    for(int j = 0; j < srcsg -> nseg; j += unitsize) {
      int tbegin = srcsg -> time[j];
      int tend = j + unitsize >= srcsg -> nseg ? srcob -> nt :
        srcsg -> time[j + unitsize];
      ret -> observset -> samples[nsample] = lrh_create_observ(srcob -> nstream,
        tend - tbegin, srcob -> ndim);
      ret -> segset -> samples[nsample] = lrh_create_seg(srcsg -> nstream, unitsize);
      lrh_observ* dstob = ret -> observset -> samples[nsample];
      lrh_seg* dstsg = ret -> segset -> samples[nsample];

      // copy segmentation
      for(int k = 0; k < unitsize; k ++) {
        dstsg -> time[k] = srcsg -> time[j + k] - tbegin;
        dstsg -> durstate[k] = srcsg -> durstate[j + k];
      }
      for(int l = 0; l < srcsg -> nstream; l ++)
        for(int k = 0; k < unitsize; k ++)
          dstsg -> outstate[l][k] = srcsg -> outstate[l][j + k];

      // copy observation
      for(int l = 0; l < srcsg -> nstream; l ++)
        for(int t = 0; t < tend - tbegin; t ++)
          for(int n = 0; n < dstob -> ndim[l]; n ++)
            lrh_obm(dstob, t, n, l) = lrh_obm(srcob, t + tbegin, n, l);

      nsample ++;
    }
  }
  return ret;
}

