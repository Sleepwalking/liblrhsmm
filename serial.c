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

#include "serial.h"

int lrh_write_duration(cmp_ctx_t* dst, lrh_duration* src) {
  if(! cmp_write_array(dst, 6)) return 0;
  if(! cmp_write_float(dst, src -> mean)) return 0;
  if(! cmp_write_float(dst, src -> var)) return 0;
  if(! cmp_write_int(dst, src -> floor)) return 0;
  if(! cmp_write_int(dst, src -> ceil)) return 0;
  if(! cmp_write_int(dst, src -> fixed_mean)) return 0;
  if(! cmp_write_float(dst, src -> vfloor)) return 0;
  return 1;
}

int lrh_write_gmm(cmp_ctx_t* dst, lrh_gmm* src) {
  if(! cmp_write_array(dst, 4)) return 0;
  if(! cmp_write_int(dst, src -> nmix)) return 0;
  if(! cmp_write_int(dst, src -> ndim)) return 0;
  
  // array for weight
  if(! cmp_write_array(dst, src -> nmix)) return 0;
  for(int i = 0; i < src -> nmix; i ++)
    if(! cmp_write_float(dst, src -> weight[i])) return 0;
  
  // array for mean & var
  if(! cmp_write_array(dst, src -> nmix * src -> ndim * 2)) return 0;
  for(int i = 0; i < src -> nmix; i ++)
    for(int j = 0; j < src -> ndim; j ++) {
      if(! cmp_write_float(dst, lrh_gmmu(src, i, j))) return 0;
      if(! cmp_write_float(dst, lrh_gmmv(src, i, j))) return 0;
    }
  
  return 1;
}

int lrh_write_stream(cmp_ctx_t* dst, lrh_stream* src) {
  if(! cmp_write_array(dst, 2)) return 0;
  if(! cmp_write_float(dst, src -> weight)) return 0;

  if(! cmp_write_array(dst, src -> ngmm)) return 0;
  for(int i = 0; i < src -> ngmm; i ++)
    if(! lrh_write_gmm(dst, src -> gmms[i])) return 0;
  return 1;
}

int lrh_write_model(cmp_ctx_t* dst, lrh_model* src) {
  if(! cmp_write_array(dst, 2)) return 0;
  
  if(! cmp_write_array(dst, src -> nstream)) return 0;
  for(int i = 0; i < src -> nstream; i ++)
    if(! lrh_write_stream(dst, src -> streams[i])) return 0;
  
  if(! cmp_write_array(dst, src -> nduration)) return 0;
  for(int i = 0; i < src -> nduration; i ++)
    if(! lrh_write_duration(dst, src -> durations[i])) return 0;
  return 1;
}

lrh_duration* lrh_read_duration(cmp_ctx_t* src) {
  uint32_t structsize;
  if(! cmp_read_array(src, & structsize)) return NULL;
  if(structsize != 6) return NULL;
  
  lrh_duration* ret = lrh_create_duration();
  float tmp;
  int32_t tmp_int;
  if(! cmp_read_float(src, & tmp)) { lrh_delete_duration(ret); return NULL;}
  ret -> mean = tmp;
  if(! cmp_read_float(src, & tmp)) { lrh_delete_duration(ret); return NULL;}
  ret -> var = tmp;
  if(! cmp_read_int(src, & tmp_int)) { lrh_delete_duration(ret); return NULL;}
  ret -> floor = tmp_int;
  if(! cmp_read_int(src, & tmp_int)) { lrh_delete_duration(ret); return NULL;}
  ret -> ceil = tmp_int;
  if(! cmp_read_int(src, & tmp_int)) { lrh_delete_duration(ret); return NULL;}
  ret -> fixed_mean = tmp_int;
  if(! cmp_read_float(src, & tmp)) { lrh_delete_duration(ret); return NULL;}
  ret -> vfloor = tmp;
  
  return ret;
}

lrh_gmm* lrh_read_gmm(cmp_ctx_t* src) {
  uint32_t structsize, arrsize;
  int32_t nmix, ndim;
  if(! cmp_read_array(src, & structsize)) return NULL;
  if(structsize != 4) return NULL;
  
  if(! cmp_read_int(src, & nmix)) return NULL;
  if(! cmp_read_int(src, & ndim)) return NULL;
  
  lrh_gmm* ret = lrh_create_gmm(nmix, ndim);
  if(! cmp_read_array(src, & arrsize)) { lrh_delete_gmm(ret); return NULL;}
  if(arrsize != nmix) { lrh_delete_gmm(ret); return NULL;}
  for(int i = 0; i < nmix; i ++) {
    float tmp;
    if(! cmp_read_float(src, & tmp)) { lrh_delete_gmm(ret); return NULL;}
    ret -> weight[i] = tmp;
  }
  
  uint32_t matsize;
  if(! cmp_read_array(src, & matsize)) { lrh_delete_gmm(ret); return NULL;}
  if(matsize != nmix * ndim * 2) { lrh_delete_gmm(ret); return NULL;}
  for(int i = 0; i < nmix; i ++)
    for(int j = 0; j < ndim; j ++) {
      float tmp;
      if(! cmp_read_float(src, & tmp)) { lrh_delete_gmm(ret); return NULL;}
      lrh_gmmu(ret, i, j) = tmp;
      if(! cmp_read_float(src, & tmp)) { lrh_delete_gmm(ret); return NULL;}
      lrh_gmmv(ret, i, j) = tmp;
    }
  
  return ret;
}

lrh_stream* lrh_read_stream(cmp_ctx_t* src) {
  uint32_t structsize, ngmm;
  if(! cmp_read_array(src, & structsize)) return NULL;
  if(structsize != 2) return NULL;
  
  float tmp;
  if(! cmp_read_float(src, & tmp)) return NULL;
  if(! cmp_read_array(src, & ngmm)) return NULL;
  lrh_stream* ret = lrh_create_empty_stream(ngmm);
  ret -> weight = tmp;
  
  for(int i = 0; i < ngmm; i ++) {
    ret -> gmms[i] = lrh_read_gmm(src);
    if(ret -> gmms[i] == NULL) { lrh_delete_stream(ret); return NULL;}
  }
  
  return ret;
}

lrh_model* lrh_read_model(cmp_ctx_t* src) {
  uint32_t structsize;
  if(! cmp_read_array(src, & structsize)) return NULL;
  if(structsize != 2) return NULL;
  
  uint32_t nstream, nduration;
  if(! cmp_read_array(src, & nstream)) return NULL;
  lrh_model* ret = lrh_create_empty_model(nstream, 1);
  
  for(int i = 0; i < nstream; i ++) {
    ret -> streams[i] = lrh_read_stream(src);
    if(ret -> streams[i] == NULL) { lrh_delete_model(ret); return NULL;}
  }
  
  if(! cmp_read_array(src, & nduration)) { lrh_delete_model(ret); return NULL;}
  ret -> nduration = nduration;
  ret -> durations = realloc(ret -> durations, nduration * sizeof(lrh_duration));
  for(int i = 0; i < nduration; i ++) {
    ret -> durations[i] = lrh_read_duration(src);
    if(ret -> durations[i] == NULL) { lrh_delete_model(ret); return NULL;}
  }
  
  return ret;
}

int lrh_write_observ(cmp_ctx_t* dst, lrh_observ* src) {
  if(! cmp_write_array(dst, 3)) return 0;
  
  if(! cmp_write_int(dst, src -> nt)) return 0;
  
  if(! cmp_write_array(dst, src -> nstream)) return 0;
  for(int i = 0; i < src -> nstream; i ++)
    if(! cmp_write_int(dst, src -> ndim[i])) return 0;
  
  if(! cmp_write_array(dst, src -> nstream)) return 0;
  for(int l = 0; l < src -> nstream; l ++) {
    if(! cmp_write_array(dst, src -> ndim[l] * src -> nt)) return 0;
    for(int t = 0; t < src -> nt; t ++)
      for(int i = 0; i < src -> ndim[l]; i ++)
        if(! cmp_write_float(dst, lrh_obm(src, t, i, l))) return 0;
  }
  
  return 1;
}

int lrh_write_observset(cmp_ctx_t* dst, lrh_observset* src) {
  if(! cmp_write_array(dst, src -> nsample)) return 0;
  for(int i = 0; i < src -> nsample; i ++)
    if(! lrh_write_observ(dst, src -> samples[i])) return 0;
  
  return 1;
}

int lrh_write_seg(cmp_ctx_t* dst, lrh_seg* src) {
  if(! cmp_write_array(dst, 3)) return 0;
  
  if(! cmp_write_array(dst, src -> nseg)) return 0;
  for(int i = 0; i < src -> nseg; i ++)
    if(! cmp_write_int(dst, src -> time[i])) return 0;
  
  if(! cmp_write_array(dst, src -> nseg)) return 0;
  for(int i = 0; i < src -> nseg; i ++)
    if(! cmp_write_int(dst, src -> durstate[i])) return 0;
  
  if(! cmp_write_array(dst, src -> nstream)) return 0;
  for(int l = 0; l < src -> nstream; l ++) {
    if(! cmp_write_array(dst, src -> nseg)) return 0;
    for(int i = 0; i < src -> nseg; i ++)
      if(! cmp_write_int(dst, src -> outstate[l][i])) return 0;
  }
  
  return 1;
}

int lrh_write_segset(cmp_ctx_t* dst, lrh_segset* src) {
  if(! cmp_write_array(dst, src -> nsample)) return 0;
  for(int i = 0; i < src -> nsample; i ++)
    if(! lrh_write_seg(dst, src -> samples[i])) return 0;
  
  return 1;
}

lrh_observ* lrh_read_observ(cmp_ctx_t* src) {
  uint32_t structsize, nstream;
  int32_t nt;
  if(! cmp_read_array(src, & structsize)) return NULL;
  if(structsize != 3) return NULL;
  
  if(! cmp_read_int(src, & nt)) return NULL;
  
  if(! cmp_read_array(src, & nstream)) return NULL;
  if(nstream <= 0) return NULL;
  int* ndim = calloc(nstream, sizeof(int));
  for(int i = 0; i < nstream; i ++)
    if(! cmp_read_int(src, & ndim[i])) { free(ndim); return NULL;}
  lrh_observ* ret = lrh_create_observ(nstream, nt, ndim);
  free(ndim);
  
  if(! cmp_read_array(src, & nstream)) { lrh_delete_observ(ret); return NULL;}
  if(nstream != ret -> nstream) { lrh_delete_observ(ret); return NULL;}
  for(int l = 0; l < nstream; l ++) {
    uint32_t lsize;
    if(! cmp_read_array(src, & lsize)) { lrh_delete_observ(ret); return NULL;}
    if(lsize != ret -> ndim[l] * nt) { lrh_delete_observ(ret); return NULL;}
    for(int t = 0; t < ret -> nt; t ++)
      for(int i = 0; i < ret -> ndim[l]; i ++) {
        float tmp;
        if(! cmp_read_float(src, & tmp)) { lrh_delete_observ(ret); return NULL;}
        lrh_obm(ret, t, i, l) = tmp;
      }
  }
  
  return ret;
}

lrh_observset* lrh_read_observset(cmp_ctx_t* src) {
  uint32_t nsample;
  if(! cmp_read_array(src, & nsample)) return NULL;
  
  lrh_observset* ret = lrh_create_empty_observset(nsample);
  for(int i = 0; i < nsample; i ++) {
    ret -> samples[i] = lrh_read_observ(src);
    if(ret -> samples[i] == NULL) { lrh_delete_observset(ret); return NULL;};
  }
  
  return ret;
}

lrh_seg* lrh_read_seg(cmp_ctx_t* src) {
  uint32_t structsize, nseg, nstream;
  if(! cmp_read_array(src, & structsize)) return NULL;
  if(structsize != 3) return NULL;
  
  if(! cmp_read_array(src, & nseg)) return NULL;
  lrh_seg* ret = malloc(sizeof(lrh_seg));
  ret -> nseg = nseg; ret -> nstream = 0; ret -> outstate = NULL;
  ret -> time = calloc(nseg, sizeof(int));
  ret -> durstate = calloc(nseg, sizeof(int));
  for(int i = 0; i < nseg; i ++)
    if(! cmp_read_int(src, & ret -> time[i])) { lrh_delete_seg(ret); return NULL;}
  
  if(! cmp_read_array(src, & nseg)) { lrh_delete_seg(ret); return NULL;}
  if(nseg != ret -> nseg) { lrh_delete_seg(ret); return NULL;}
  for(int i = 0; i < nseg; i ++)
    if(! cmp_read_int(src, & ret -> durstate[i])) { lrh_delete_seg(ret); return NULL;}
  
  if(! cmp_read_array(src, & nstream)) { lrh_delete_seg(ret); return NULL;}
  ret -> nstream = nstream;
  ret -> outstate = calloc(nstream, sizeof(int*));
  for(int l = 0; l < nstream; l ++) {
    ret -> outstate[l] = calloc(nseg, sizeof(int));
    if(! cmp_read_array(src, & nseg)) { lrh_delete_seg(ret); return NULL;}
    for(int i = 0; i < nseg; i ++)
      if(! cmp_read_int(src, & ret -> outstate[l][i])) { lrh_delete_seg(ret); return NULL;}
  }
  
  return ret;
}

lrh_segset* lrh_read_segset(cmp_ctx_t* src) {
  uint32_t nsample;
  if(! cmp_read_array(src, & nsample)) return NULL;
  
  lrh_segset* ret = lrh_create_empty_segset(nsample);
  for(int i = 0; i < nsample; i ++) {
    ret -> samples[i] = lrh_read_seg(src);
    if(ret -> samples[i] == NULL) { lrh_delete_segset(ret); return NULL;};
  }
  
  return ret;
}

