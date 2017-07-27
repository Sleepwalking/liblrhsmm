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

#include <stdio.h>
#include "../estimate.h"
#include "../generate.h"
#include "../inference.h"
#include "../math-funcs.h"
#include "../mempool.h"
#include "../serial.h"
#include "test-common.h"

static void inference_sample(lrh_model* h, lrh_observ* ob, lrh_seg* seg) {
  lrh_mempool* pool = lrh_create_mempool(1024 * 1024);
  FP_TYPE* outp   = lrh_sample_outputprob_lg(h, ob, seg);
  int* reseg      = lrh_viterbi (h, seg, outp, ob -> nt, NULL);
  FP_TYPE* ab, *bb;
  FP_TYPE* a      = lrh_forward (h, seg, outp, ob -> nt, & ab);
  FP_TYPE* av     = lrh_forward_geometric(h, seg, outp, ob -> nt);
  FP_TYPE* b      = lrh_backward(h, seg, outp, ob -> nt, & bb);
  FP_TYPE* bv     = lrh_backward_geometric(h, seg, outp, ob -> nt);
  FP_TYPE total   = lrh_total_backward_b(bb, outp, seg -> nseg);
  FP_TYPE totalv  = lrh_total_pseudo(av, ob -> nt, seg -> nseg);
  lrh_pslice* dp  = lrh_durocp  (h, seg, a, b, outp, total, ob -> nt, pool);
  FP_TYPE* sp     = lrh_stocp   (h, seg, dp, ob -> nt);
  FP_TYPE* spv    = lrh_stocp_geometric(h, seg, av, bv, totalv, ob -> nt);
  lrh_pslice* tpv = lrh_sttran_geometric(h, seg, av, bv, outp, totalv, ob -> nt, pool);

  /*
  for(int t = 0; t < ob -> nt; t ++) {
    for(int i = 0; i < seg -> nseg; i ++) {
      if(tpv[t * seg -> nseg + i].upper > 1)
        fprintf(stderr, "%f,", sp[t * seg -> nseg + i]);
      else
        fprintf(stderr, "0,");
    }
    fprintf(stderr, "\n");
  }
  exit(1);*/

  printf("Re-alignment:\n");
  for(int i = 0; i < seg -> nseg; i ++)
    printf("%d\t%d\t(%d)\n", reseg[i * 2 + 0], seg -> time[i],
      reseg[i * 2 + 0] - seg -> time[i]);

  printf("Total probabilities estimated at each state:\n");
  for(int i = 0; i < seg -> nseg; i ++) {
    FP_TYPE total_i = lrh_total(a, b, ob -> nt, seg -> nseg, i);
    printf("total (%d) = %f\n", i, total_i);
    if(fabs(total_i - total) > fabs(total / 1e5)) {
      fprintf(stderr, "Warning: total probability estimated from both forward and backward probability "
        "is different from the one estimated from backward probability only (and the difference is beyond "
        "the typical range of numerical error).\n");
      //exit(1);
    }
  }
  
  printf("total (backward) = %f\n", total);
  printf("total (pseudo) = %f\n", lrh_total_pseudo(a, ob -> nt, seg -> nseg));
  printf("total (geometric forward) = %f\n", lrh_total_pseudo(av, ob -> nt, seg -> nseg));
  printf("total (geometric backward) = %f\n", lrh_total_pseudo_backward(bv, outp, seg -> nseg));

  free(a); free(b); free(ab); free(bb);
  free(dp); free(sp); free(av); free(bv); free(spv); free(tpv); free(outp);
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
    lrh_seg_buildjumps(dataset_seg -> samples[i]);
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

