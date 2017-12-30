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

#include <stdlib.h>
#include <string.h>
#include "math-funcs.h"
#include "model.h"
#include "common.h"

#define gmmu lrh_gmmu
#define gmmv lrh_gmmv

lrh_gmm* lrh_create_gmm(int nmix, int ndim) {
  lrh_gmm* ret = malloc(sizeof(lrh_gmm));
  ret -> mean = calloc(nmix * ndim, sizeof(FP_TYPE));
  ret -> var = calloc(nmix * ndim, sizeof(FP_TYPE));
  ret -> weight = calloc(nmix, sizeof(FP_TYPE));

  ret -> vfloor = calloc(nmix * ndim, sizeof(FP_TYPE));
  ret -> _tmp_term = calloc(nmix, sizeof(FP_TYPE));
  ret -> nmix = nmix;
  ret -> ndim = ndim;
  for(int i = 0; i < nmix; i ++) {
    ret -> weight[i] = 1.0 / nmix;
    ret -> _tmp_term[i] = LRH_UNINITIALIZED;
  }
  return ret;
}

lrh_duration* lrh_create_duration() {
  lrh_duration* ret = malloc(sizeof(lrh_duration));
  ret -> mean = 0;
  ret -> var = 0;
  ret -> floor = -1;
  ret -> ceil = -1;
  ret -> fixed_mean = -1;
  ret -> vfloor = 0;
  ret -> _tmp_term = LRH_UNINITIALIZED;
  ret -> _tmp_prep = NULL;
  return ret;
}

lrh_stream* lrh_create_empty_stream(int ngmm) {
  lrh_stream* ret = malloc(sizeof(lrh_stream));
  ret -> gmms = calloc(ngmm, sizeof(lrh_gmm*));
  ret -> ngmm = ngmm;
  ret -> weight = 1;
  return ret;
}

lrh_stream* lrh_create_stream(int ngmm, int nmix, int ndim) {
  lrh_stream* ret = lrh_create_empty_stream(ngmm);
  for(int i = 0; i < ngmm; i ++)
    ret -> gmms[i] = lrh_create_gmm(nmix, ndim);
  return ret;
}

lrh_model* lrh_create_empty_model(int nstream, int nduration) {
  lrh_model* ret = malloc(sizeof(lrh_model));
  ret -> nstream = nstream;
  ret -> nduration = nduration;
  ret -> streams = calloc(nstream, sizeof(lrh_stream*));
  ret -> durations = calloc(nduration, sizeof(lrh_duration*));
  return ret;
}

lrh_model* lrh_create_model(int nstream, int* ngmm, int nduration, int* nmix, int* ndim) {
  lrh_model* ret = lrh_create_empty_model(nstream, nduration);
  for(int l = 0; l < nstream; l ++) {
    ret -> streams[l] = lrh_create_empty_stream(ngmm[l]);
    for(int i = 0; i < ngmm[l]; i ++)
      ret -> streams[l] -> gmms[i] = lrh_create_gmm(nmix[l], ndim[l]);
  }
  for(int i = 0; i < nduration; i ++)
    ret -> durations[i] = lrh_create_duration();
  return ret;
}

void lrh_delete_gmm(lrh_gmm* dst) {
  if(dst == NULL) return;
  free(dst -> mean);
  free(dst -> var);
  free(dst -> vfloor);
  free(dst -> weight);
  free(dst -> _tmp_term);
  free(dst);
}

void lrh_delete_duration(lrh_duration* dst) {
  if(dst == NULL) return;
  if(dst -> _tmp_prep != NULL)
    free(dst -> _tmp_prep);
  free(dst);
}

void lrh_delete_stream(lrh_stream* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> ngmm; i ++)
    lrh_delete_gmm(dst -> gmms[i]);
  free(dst -> gmms);
  free(dst);
}

void lrh_delete_model(lrh_model* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> nstream; i ++)
    lrh_delete_stream(dst -> streams[i]);
  for(int i = 0; i < dst -> nduration; i ++)
    lrh_delete_duration(dst -> durations[i]);
  free(dst -> streams);
  free(dst -> durations);
  free(dst);
}

lrh_gmm* lrh_gmm_copy(lrh_gmm* src) {
  lrh_gmm* ret = lrh_create_gmm(src -> nmix, src -> ndim);
  memcpy(ret -> mean, src -> mean, sizeof(FP_TYPE) * src -> nmix * src -> ndim);
  memcpy(ret -> var , src -> var , sizeof(FP_TYPE) * src -> nmix * src -> ndim);
  memcpy(ret -> vfloor, src -> vfloor, sizeof(FP_TYPE) * src -> nmix * src -> ndim);
  memcpy(ret -> weight, src -> weight, sizeof(FP_TYPE) * src -> nmix);
  return ret;
}

FP_TYPE lrh_duration_logp(lrh_duration* src, FP_TYPE duration) {
  if(src -> _tmp_term == LRH_UNINITIALIZED)
    src -> _tmp_term = -log(src -> var * 2.0 * M_PI) * 0.5;
  if(src -> floor > 0 && duration < src -> floor) return NEGINF;
  FP_TYPE p = src -> _tmp_term - 
    (duration - src -> mean) * (duration - src -> mean) / src -> var / 2.0;
  return p * lrh_inference_duration_weight;
}

lrh_gmm* lrh_combine_duration(lrh_duration** src, int nsrc) {
  lrh_gmm* ret = lrh_create_gmm(1, nsrc);
  for(int i = 0; i < nsrc; i ++) {
    ret -> mean[i] = src[i] -> mean;
    ret -> var[i]  = src[i] -> var;
  }
  return ret;
}

lrh_gmm* lrh_combine_gmm(lrh_gmm** src, int nsrc) {
  int ndim = 0;
  for(int i = 0; i < nsrc; i ++) {
    ndim += src[i] -> ndim;
    if(src[i] -> nmix != 1)
      return NULL; // doesn't support multiple-mixture gmm
  }
  lrh_gmm* ret = lrh_create_gmm(1, ndim);
  int n = 0;
  for(int i = 0; i < nsrc; i ++)
    for(int j = 0; j < src[i] -> ndim; j ++) {
      ret -> mean[n] = src[i] -> mean[j];
      ret -> vfloor[n] = src[i] -> vfloor[j];
      ret -> var[n] = src[i] -> var[j];
      n ++;
    }
  return ret;
}

static void lrh_gmm_mixuniinc(lrh_gmm* dst, FP_TYPE d) {
  int nsrc = 0;
  FP_TYPE maxw = 0;
  for(int i = 0; i < dst -> nmix; i ++)
    if(dst -> weight[i] > maxw) {
      maxw = dst -> weight[i];
      nsrc = i;
    }
  int ndst = dst -> nmix;
  dst -> weight[ndst] = dst -> weight[nsrc] * 0.5;
  dst -> weight[nsrc] *= 0.5;
  for(int i = 0; i < dst -> ndim; i ++) {
    lrh_gmmu(dst, ndst, i) = lrh_gmmu(dst, nsrc, i) + sqrt(lrh_gmmv(dst, nsrc, i)) * d;
    lrh_gmmv(dst, ndst, i) = lrh_gmmv(dst, nsrc, i);
    lrh_gmmu(dst, nsrc, i) -= sqrt(lrh_gmmv(dst, nsrc, i)) * d;
    lrh_gmmvf(dst, ndst, i) = lrh_gmmvf(dst, nsrc, i);
  }
  dst -> nmix ++;
}

void lrh_gmm_mixinc(lrh_gmm* dst, int nmix, FP_TYPE d) {
  if(dst -> nmix >= nmix) return;
  dst -> mean = realloc(dst -> mean, sizeof(FP_TYPE) * nmix * dst -> ndim);
  dst -> var  = realloc(dst -> var , sizeof(FP_TYPE) * nmix * dst -> ndim);
  dst -> vfloor = realloc(dst -> vfloor, sizeof(FP_TYPE) * nmix * dst -> ndim);
  dst -> weight = realloc(dst -> weight, sizeof(FP_TYPE) * nmix);
  while(dst -> nmix != nmix)
    lrh_gmm_mixuniinc(dst, d);
/*
  for(int i = 0; i < dst -> nmix; i ++) {
    if(isnan(dst -> weight[i]))
      printf("%d mixture, weight is nan.\n", i);
    for(int j = 0; j < dst -> ndim; j ++)
      if(isnan(lrh_gmmu(dst, i, j)))
        printf("%d mixture, mean[%d] is nan.\n", i, j);
    for(int j = 0; j < dst -> ndim; j ++)
      if(isnan(lrh_gmmv(dst, i, j)))
        printf("%d mixture, var[%d] is nan.\n", i, j);
  }
*/
}

void lrh_model_precompute(lrh_model* dst) {
# pragma omp critical
  for(int l = 0; l < dst -> nstream; l ++) {
    for(int i = 0; i < dst -> streams[l] -> ngmm; i ++) {
      lrh_gmm* igmm = dst -> streams[l] -> gmms[i];
      if(igmm -> nmix == 0) continue;
      if(igmm -> _tmp_term[0] == LRH_UNINITIALIZED) {
        for(int k = 0; k < igmm -> nmix; k ++) {
          FP_TYPE mul = 0;
          for(int j = 0; j < igmm -> ndim; j ++)
            mul -= log(lrh_gmmv(igmm, k, j) * 2.0 * M_PI);
          igmm -> _tmp_term[k] = mul;
        }
      }
    }
  }
# pragma omp critical
  for(int i = 0; i < dst -> nduration; i ++) {
    lrh_duration* idur = dst -> durations[i];
    if(idur -> _tmp_term == LRH_UNINITIALIZED) {
      idur -> _tmp_prep = realloc(idur -> _tmp_prep, lrh_precompute_duration * sizeof(FP_TYPE));
      for(int j = 0; j < lrh_precompute_duration; j ++)
        idur -> _tmp_prep[j] = lrh_duration_logp(idur, j) * lrh_daem_temperature;
    }
  }
}

lrh_gmm_stat* lrh_gmm_stat_from_gmm(lrh_gmm* src) {
  lrh_gmm_stat* ret = malloc(sizeof(lrh_gmm_stat));
  ret -> valsum = calloc(src -> nmix * src -> ndim, sizeof(FP_TYPE));
  ret -> sqrsum = calloc(src -> nmix * src -> ndim, sizeof(FP_TYPE));
  ret -> weightsum = calloc(src -> nmix, sizeof(FP_TYPE));
  ret -> nmix = src -> nmix;
  ret -> ndim = src -> ndim;
  lrh_gmm_stat_clear(ret);
  return ret;
}

lrh_duration_stat* lrh_duration_stat_from_duration(lrh_duration* src) {
  lrh_duration_stat* ret = malloc(sizeof(lrh_duration_stat));
  lrh_duration_stat_clear(ret);
  return ret;
}

lrh_stream_stat* lrh_stream_stat_from_stream(lrh_stream* src) {
  lrh_stream_stat* ret = malloc(sizeof(lrh_stream_stat));
  ret -> gmms = calloc(src -> ngmm, sizeof(lrh_gmm_stat*));
  ret -> ngmm = src -> ngmm;
  for(int i = 0; i < src -> ngmm; i ++)
    ret -> gmms[i] = lrh_gmm_stat_from_gmm(src -> gmms[i]);
  return ret;
}

lrh_model_stat* lrh_model_stat_from_model(lrh_model* src) {
  lrh_model_stat* ret = malloc(sizeof(lrh_model_stat));
  ret -> streams = calloc(src -> nstream, sizeof(lrh_stream_stat*));
  ret -> durations = calloc(src -> nduration, sizeof(lrh_duration_stat*));
  ret -> nstream = src -> nstream;
  ret -> nduration = src -> nduration;
  for(int i = 0; i < src -> nstream; i ++)
    ret -> streams[i] = lrh_stream_stat_from_stream(src -> streams[i]);
  for(int i = 0; i < src -> nduration; i ++)
    ret -> durations[i] = lrh_duration_stat_from_duration(src -> durations[i]);
  return ret;
}

void lrh_delete_gmm_stat(lrh_gmm_stat* dst) {
  if(dst == NULL) return;
  free(dst -> valsum);
  free(dst -> sqrsum);
  free(dst -> weightsum);
  free(dst);
}

void lrh_delete_duration_stat(lrh_duration_stat* dst) {
  if(dst == NULL) return;
  free(dst);
}

void lrh_delete_stream_stat(lrh_stream_stat* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> ngmm; i ++)
    lrh_delete_gmm_stat(dst -> gmms[i]);
  free(dst -> gmms);
  free(dst);
}

void lrh_delete_model_stat(lrh_model_stat* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> nduration; i ++)
    lrh_delete_duration_stat(dst -> durations[i]);
  free(dst -> durations);
  for(int i = 0; i < dst -> nstream; i ++)
    lrh_delete_stream_stat(dst -> streams[i]);
  free(dst -> streams);
  free(dst);
}

void lrh_gmm_stat_clear(lrh_gmm_stat* dst) {
  memset(dst -> valsum, 0, sizeof(FP_TYPE) * dst -> nmix * dst -> ndim);
  memset(dst -> sqrsum, 0, sizeof(FP_TYPE) * dst -> nmix * dst -> ndim);
  memset(dst -> weightsum, 0, sizeof(FP_TYPE) * dst -> nmix);
  dst -> ocpsum = 0;
}

void lrh_duration_stat_clear(lrh_duration_stat* dst) {
  dst -> valsum = 0;
  dst -> sqrsum = 0;
  dst -> ocpsum = 0;
}

void lrh_stream_stat_clear(lrh_stream_stat* dst) {
  for(int i = 0; i < dst -> ngmm; i ++)
    lrh_gmm_stat_clear(dst -> gmms[i]);
}

void lrh_model_stat_clear(lrh_model_stat* dst) {
  for(int i = 0; i < dst -> nstream; i ++)
    lrh_stream_stat_clear(dst -> streams[i]);
  for(int i = 0; i < dst -> nduration; i ++)
    lrh_duration_stat_clear(dst -> durations[i]);
}
