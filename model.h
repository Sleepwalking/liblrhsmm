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

#ifndef LRHSMM_MODEL_H
#define LRHSMM_MODEL_H

typedef struct {
  FP_TYPE mean, var; // mean and variance
  int floor;         // minimum duration
  int ceil;          // maximum duration
  int fixed_mean;
  FP_TYPE vfloor;

  FP_TYPE _tmp_term; // temporary variable for storing the constant term in normal pdf
  FP_TYPE* _tmp_prep;  // pre-calculated duration probabilities (whose length is lrh_precompute_duration)
} lrh_duration;

typedef struct {
  FP_TYPE valsum, sqrsum; // sum of values and sum of squared values
  FP_TYPE ocpsum;         // sum of occupancy
} lrh_duration_stat;

typedef struct {
  FP_TYPE* mean;   // matrix of the mean vector in each mixture
  FP_TYPE* var;    // matrix of the diagonal variance vector in each mixture
  FP_TYPE* weight; // vector of mixture weights
  FP_TYPE* vfloor; // matrix of variance floor, associated with var

  FP_TYPE* _tmp_term; // vector of constant terms (one for each mixture component)
  int nmix, ndim;  // number of mixture and dimension
} lrh_gmm;

typedef struct {
  FP_TYPE* valsum;    // matrix of the sum of weighted observation vectors in each mixture
  FP_TYPE* sqrsum;    // matrix of the sum of weighted element-wise-squared observation vectors in each mixture
  FP_TYPE* weightsum; // vector of the sum of weighted mixture occupancy
  FP_TYPE ocpsum;     // sum of state occupancy

  int nmix, ndim;
} lrh_gmm_stat;

#define lrh_gmmu(obj, k, i) (((obj) -> mean)[k * ((obj) -> ndim) + i]) // access mean of mixture k, dimension i
#define lrh_gmmv(obj, k, i) (((obj) -> var )[k * ((obj) -> ndim) + i]) // access variance of mixture k, dimension i
#define lrh_gmmvf(obj, k, i) (((obj) -> vfloor)[k * ((obj) -> ndim) + i]) // access variance floor
#define lrh_gmmvs(obj, k, i) (((obj) -> valsum)[k * ((obj) -> ndim) + i]) // access 1st momentum sum
#define lrh_gmmss(obj, k, i) (((obj) -> sqrsum)[k * ((obj) -> ndim) + i]) // access 2nd momentum sum

typedef struct {
  lrh_gmm** gmms; // array of Gaussian mixture models
  int ngmm;       // number of Gaussian mixture models
  FP_TYPE weight;  // weight of this stream
} lrh_stream;

typedef struct {
  lrh_gmm_stat** gmms;
  int ngmm;
} lrh_stream_stat;

typedef struct {
  lrh_stream** streams;     // array of streams (each containing a group of gmms)
  lrh_duration** durations; // array of duration pdfs
  int nstream, nduration;
} lrh_model;

typedef struct {
  lrh_stream_stat** streams;
  lrh_duration_stat** durations;
  int nstream, nduration;
} lrh_model_stat;

lrh_gmm* lrh_create_gmm(int nmix, int ndim);
lrh_duration* lrh_create_duration();
lrh_stream* lrh_create_stream(int ngmm, int nmix, int ndim);
lrh_stream* lrh_create_empty_stream(int ngmm);
lrh_model* lrh_create_model(int nstream, int* ngmm, int nduration, int* nmix, int* ndim);
lrh_model* lrh_create_empty_model(int nstream, int nduration);
void lrh_delete_gmm(lrh_gmm* dst);
void lrh_delete_duration(lrh_duration* dst);
void lrh_delete_stream(lrh_stream* dst);
void lrh_delete_model(lrh_model* dst);
lrh_gmm* lrh_gmm_copy(lrh_gmm* src);
FP_TYPE lrh_duration_logp(lrh_duration* src, FP_TYPE duration); // without using pre-computed value

// combine multiple lrh_durations into one single-mixture gmm with dimension nsrc
lrh_gmm* lrh_combine_duration(lrh_duration** src, int nsrc);
// combine multiple single-mixture lrh_gmms into one single-mixture gmm with dimensions stacked together
lrh_gmm* lrh_combine_gmm(lrh_gmm** src, int nsrc);
// increase the number of mixture components by iteratively splitting the largest component
void lrh_gmm_mixinc(lrh_gmm* dst, int nmix, FP_TYPE d);
// calculate constant values before doing inference
void lrh_gmm_precompute(lrh_gmm* dst);
void lrh_model_precompute(lrh_model* dst);

lrh_gmm_stat* lrh_gmm_stat_from_gmm(lrh_gmm* src);
lrh_duration_stat* lrh_duration_stat_from_duration(lrh_duration* src);
lrh_stream_stat* lrh_stream_stat_from_stream(lrh_stream* src);
lrh_model_stat* lrh_model_stat_from_model(lrh_model* src);
void lrh_delete_gmm_stat(lrh_gmm_stat* dst);
void lrh_delete_duration_stat(lrh_duration_stat* dst);
void lrh_delete_stream_stat(lrh_stream_stat* dst);
void lrh_delete_model_stat(lrh_model_stat* dst);

// reset expectation collectors
void lrh_gmm_stat_clear(lrh_gmm_stat* dst);
void lrh_duration_stat_clear(lrh_duration_stat* dst);
void lrh_stream_stat_clear(lrh_stream_stat* dst);
void lrh_model_stat_clear(lrh_model_stat* dst);

#endif
