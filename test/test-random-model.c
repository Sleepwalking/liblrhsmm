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

#include <stdio.h>
#include "../estimate.h"
#include "../generate.h"
#include "../inference.h"
#include "../math-funcs.h"
#include "../mempool.h"
#include "../serial.h"

inline static bool read_bytes(void* data, size_t sz, FILE* fh) {
  return fread(data, sizeof(uint8_t), sz, fh) == (sz * sizeof(uint8_t));
}

inline static bool file_reader(cmp_ctx_t* ctx, void* data, size_t limit) {
  return read_bytes(data, limit, (FILE*)ctx -> buf);
}

inline static size_t file_writer(cmp_ctx_t* ctx, const void* data, size_t count) {
  return fwrite(data, sizeof(uint8_t), count, (FILE*)ctx -> buf);
}


static lrh_model* gen_random_model(int nstream, int nstate, int nmix, int dims[]) {
  lrh_model* h = lrh_create_empty_model(nstream, nstate);
  h -> streams[0] = lrh_create_stream(nstate, nmix, dims[0]);
  h -> streams[1] = lrh_create_stream(nstate, nmix, dims[1]);
  
  for(int l = 0; l < nstream; l ++)
    for(int i = 0; i < nstate; i ++)
      for(int k = 0; k < nmix; k ++)
        for(int j = 0; j < h -> streams[l] -> gmms[i] -> ndim; j ++) {
          lrh_gmmu(h -> streams[l] -> gmms[i], k, j) = lrh_random() * 5 - 2.5;
          lrh_gmmv(h -> streams[l] -> gmms[i], k, j) = lrh_random() * 1 + 3;
        }
  for(int i = 0; i < nstate; i ++) {
    h -> durations[i] = lrh_create_duration();
    h -> durations[i] -> mean = lrh_random() * 15 + 5;
    h -> durations[i] -> var = lrh_random() * 2 + 3;
  }
  return h;
}

static lrh_seg* gen_random_seq(int nstream, int noutstates[], int ndurstate, int nseg) {
  lrh_seg* ret = lrh_create_seg(nstream, nseg);
  for(int i = 0; i < nstream; i ++)
    for(int j = 0; j < nseg; j ++)
      ret -> outstate[i][j] = rand() % noutstates[i];
  for(int j = 0; j < nseg; j ++)
    ret -> durstate[j] = rand() % ndurstate;
  return ret;
}

static void print_model(lrh_model* h) {
  for(int i = 0; i < h -> nduration; i ++) {
    printf("duration %d/%d:\n", i, h -> nduration);
    printf("\tmean = %6.3f, var = %6.3f\n", h -> durations[i] -> mean, h -> durations[i] -> var);
  }
  for(int l = 0; l < h -> nstream; l ++) {
    printf("stream %d/%d:\n", l, h -> nstream);
    printf("\tweight = %f\n", h -> streams[l] -> weight);
    for(int i = 0; i < h -> streams[l] -> ngmm; i ++) {
      lrh_gmm* g = h -> streams[l] -> gmms[i];
      printf("\tgmm %d/%d:\n", i, h -> streams[l] -> ngmm);
      printf("\t\tnmix = %d, ndim = %d\n", g -> nmix, g -> ndim);
      for(int k = 0; k < g -> nmix; k ++) {
        printf("\t\tmixture %d/%d (weight = %f):\n", k, g -> nmix, g -> weight[k]);
        printf("\t\t\tmean:    ");
        for(int j = 0; j < g -> ndim; j ++)
          printf(" %6.3f", lrh_gmmu(g, k, j));
        printf("\n\t\t\tvariance:");
        for(int j = 0; j < g -> ndim; j ++)
          printf(" %6.3f", lrh_gmmv(g, k, j));
        printf("\n");
      }
    }
  }
}

static void print_observ(lrh_observ* ob) {
  printf("nstream = %d, nt = %d\n", ob -> nstream, ob -> nt);
  printf("dimensions:");
  for(int i = 0; i < ob -> nstream; i ++)
    printf(" %d", ob -> ndim[i]);
  printf("\n");
  for(int t = 0; t < ob -> nt; t ++) {
    printf("t = %d/%d:\n", t, ob -> nt);
    for(int l = 0; l < ob -> nstream; l ++) {
      printf("\tstream %d\n\t", l);
      for(int i = 0; i < ob -> ndim[l]; i ++)
        printf(" %f", lrh_obm(ob, t, i, l));
      printf("\n");
    }
  }
}

static void print_observ_csv(lrh_observ* ob) {
  for(int t = 0; t < ob -> nt; t ++) {
    for(int l = 0; l < ob -> nstream; l ++) {
      for(int i = 0; i < ob -> ndim[l]; i ++)
        printf("%f,", lrh_obm(ob, t, i, l));
    }
    printf("\n");
  }
}

static void inference_sample(lrh_model* h, lrh_observ* ob, lrh_seg* seg) {
  lrh_mempool* pool = lrh_create_mempool(1024 * 1024);
  FP_TYPE* outp   = lrh_sample_outputprob_lg(h, ob, seg);
  int* realign    = lrh_viterbi (h, seg, outp, ob -> nt, NULL);
  FP_TYPE* a      = lrh_forward (h, seg, outp, ob -> nt);
  FP_TYPE* av     = lrh_forward_geometric(h, seg, outp, ob -> nt);
  FP_TYPE* b      = lrh_backward(h, seg, outp, ob -> nt);
  FP_TYPE* bv     = lrh_backward_geometric(h, seg, outp, ob -> nt);
  FP_TYPE total   = lrh_total_forward(a, outp, ob -> nt, seg -> nseg);
  FP_TYPE totalv  = lrh_total_pseudo(av, ob -> nt, seg -> nseg);
  lrh_pslice* dp  = lrh_durocp  (h, seg, a, b, outp, total, ob -> nt, pool);
  FP_TYPE* sp     = lrh_stocp   (h, seg, dp, ob -> nt);
  FP_TYPE* spv    = lrh_stocp_geometric(h, seg, av, bv, totalv, ob -> nt);
  lrh_pslice* tpv = lrh_sttran_geometric(h, seg, av, bv, outp, totalv, ob -> nt, pool);
/*
  for(int t = 0; t < ob -> nt; t ++) {
    for(int i = 0; i < seg -> nseg; i ++) {
      if(tpv[t * seg -> nseg + i].upper > 1)
        fprintf(stderr, "%f,", tpv[t * seg -> nseg + i].p[1]);
      else
        fprintf(stderr, "0,");
    }
    fprintf(stderr, "\n");
  }*/

  printf("Re-alignment:\n");
  for(int i = 0; i < seg -> nseg; i ++)
    printf("%d\t%d\t(%d)\n", realign[i], seg -> time[i], realign[i] - seg -> time[i]);

  printf("Total probabilities estimated at each state:\n");
  for(int i = 0; i < seg -> nseg; i ++) {
    FP_TYPE total_i = lrh_total(a, b, ob -> nt, seg -> nseg, i);
    printf("total (%d) = %f\n", i, total_i);
    if(fabs(total_i - total) > fabs(total / 1e5)) {
      fprintf(stderr, "Error: total probability estimated from both forward and backward probability "
        "is different from the one estimated from forward probability only (and the difference is beyond "
        "the typical range of numerical error).\n");
      exit(1);
    }
  }
  
  FP_TYPE total_bk = lrh_total_backward(b, outp, seg -> nseg);
  printf("total (forward) = %f\n", total);
  printf("total (backward) = %f\n", total_bk);
  printf("total (pseudo) = %f\n", lrh_total_pseudo(a, ob -> nt, seg -> nseg));
  printf("total (geometric forward) = %f\n", lrh_total_pseudo(av, ob -> nt, seg -> nseg));
  printf("total (geometric backward) = %f\n", lrh_total_pseudo_backward(bv, outp, seg -> nseg));
  if(fabs(total_bk - total) > fabs(total / 1e5)) {
    fprintf(stderr, "Error: total probability estimated from backward probability "
      "is different from the one estimated from forward probability (and the difference is beyond "
      "the typical range of numerical error).\n");
    exit(1);
  }

  free(a); free(b); free(dp); free(sp); free(av); free(bv); free(spv); free(tpv); free(outp);
  free(realign);
  lrh_delete_mempool(pool);
}

int main() {
  int nstream = 2;
  int dims[2] = {3, 5};
  int nstate = 3;
  int nmixture = 2;
  int noutstates[2] = {3, 3};
  
  printf("Creating random model with %d stream(s), %d duration state(s) and %d mixture(s) in each state.\n",
    nstream, nstate, nmixture);
  lrh_model* h = gen_random_model(nstream, nstate, nmixture, dims);
  printf("Initial model:\n");
  print_model(h);
  
  printf("Generating random observation sequences and segmentations from the model.\n");
  int nsample = 50;
  lrh_segset* dataset_seg = lrh_create_empty_segset(nsample);
  lrh_observset* dataset_ob = lrh_create_empty_observset(nsample);
  for(int i = 0; i < nsample; i ++) {
    dataset_seg -> samples[i] = gen_random_seq(nstream, noutstates, nstate, 50);
    dataset_ob -> samples[i] = lrh_generate(h, dataset_seg -> samples[i]);
  }
  //print_observ_csv(ob);

  cmp_ctx_t cmpobj;

  printf("Writing data...\n");
  FILE* fout = fopen("rand-50.lrho", "wb");
  cmp_init(& cmpobj, fout, file_reader, file_writer);
  lrh_write_observset(& cmpobj, dataset_ob);
  fclose(fout);

  fout = fopen("rand-50.lrhs", "wb");
  cmp_init(& cmpobj, fout, file_reader, file_writer);
  lrh_write_segset(& cmpobj, dataset_seg);
  fclose(fout);
  
  lrh_delete_observset(dataset_ob);
  lrh_delete_segset(dataset_seg);
  
  printf("Loading data...\n");
  FILE* fin = fopen("rand-50.lrho", "rb");
  cmp_init(& cmpobj, fin, file_reader, file_writer);
  dataset_ob = lrh_read_observset(& cmpobj);
  fclose(fin);
  
  fin = fopen("rand-50.lrhs", "rb");
  cmp_init(& cmpobj, fin, file_reader, file_writer);
  dataset_seg = lrh_read_segset(& cmpobj);
  fclose(fin);

  printf("Dropping one of the mixtures for each GMM.\n");
  for(int l = 0; l < h -> nstream; l ++)
    for(int i = 0; i < h -> streams[l] -> ngmm; i ++) {
      lrh_gmm* dstgmm = h -> streams[l] -> gmms[i];
      dstgmm -> nmix = 1;
      dstgmm -> weight[0] = 1.0;
    }
  
  printf("Re-initializing model.\n");
  lrh_model_stat* hstat = lrh_model_stat_from_model(h);
  for(int i = 0; i < nsample; i ++)
    lrh_collect_init(hstat, dataset_ob -> samples[i], dataset_seg -> samples[i]);
  lrh_model_update(h, hstat, 0);
  lrh_delete_model_stat(hstat);

  printf("Testing inference algorithms on the initialized model.\n");
  inference_sample(h, dataset_ob -> samples[0], dataset_seg -> samples[0]);
  
  printf("Running Baum-Welch training.\n");
  FP_TYPE total_best = -1e20;
  for(int iter = 0; iter < 10; iter ++) {
    lrh_model_stat* hstat = lrh_model_stat_from_model(h);
    FP_TYPE total_lh = 0;
    for(int i = 0; i < nsample; i ++)
      total_lh += lrh_estimate(hstat, h, dataset_ob -> samples[i], dataset_seg -> samples[i]);
    if(total_lh >= total_best)
      total_best = total_lh;
    else if(total_lh < total_best - 1.0) { // taking numerical error into account
      fprintf(stderr, "Error: training doesn't converge!\n");
      exit(1);
    }
    lrh_model_update(h, hstat, 0);
    lrh_delete_model_stat(hstat);
    printf("Iteration %d, total log likelihood = %f\n", iter, total_lh);
    //return 0;
  }

  printf("Re-estimated model:\n");
  print_model(h);
  
  printf("Raising number of mixtures to 2.\n");
  for(int l = 0; l < h -> nstream; l ++)
    for(int i = 0; i < h -> streams[l] -> ngmm; i ++)
      lrh_gmm_mixinc(h -> streams[l] -> gmms[i], 2, 0.1);
  
  printf("Writing model...\n");
  fout = fopen("mix2.lrh", "wb");
  cmp_init(& cmpobj, fout, file_reader, file_writer);
  lrh_write_model(& cmpobj, h);
  fclose(fout);

  lrh_delete_model(h);
  
  printf("Loading model...\n");
  fin = fopen("mix2.lrh", "rb");
  cmp_init(& cmpobj, fin, file_reader, file_writer);
  h = lrh_read_model(& cmpobj);
  fclose(fin);

  printf("Running Baum-Welch training (2 mixtures).\n");
  total_best = -1e100;
  for(int iter = 0; iter < 20; iter ++) {
    lrh_model_stat* hstat = lrh_model_stat_from_model(h);
    FP_TYPE total_lh = 0;
    for(int i = 0; i < nsample; i ++)
      total_lh += lrh_estimate(hstat, h, dataset_ob -> samples[i], dataset_seg -> samples[i]);
    if(total_lh >= total_best)
      total_best = total_lh;
    else if(total_lh < total_best - 1.0) { // taking numerical error into account
      fprintf(stderr, "Error: training doesn't converge!\n");
      // sometimes it diverges at the first few iterations due to approximated
      //   log/exp functions
      // exit(1);
    }
    lrh_model_update(h, hstat, 0);
    lrh_delete_model_stat(hstat);
    printf("Iteration %d, total log likelihood = %f\n", iter, total_lh);
  }

  printf("Re-estimated model (2 mixtures):\n");
  print_model(h);

  lrh_delete_model(h);
  lrh_delete_segset(dataset_seg);
  lrh_delete_observset(dataset_ob);
  return 0;
}

