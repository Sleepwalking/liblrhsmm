#pragma once

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

static void print_seg(lrh_seg* seg) {
  printf("nseg = %d\n", seg -> nseg);
  for(int i = 0; i < seg -> nseg; i ++) {
    printf("\tt = %d\n", seg -> time[i]);
    printf("\tdurstate = %d\n", seg -> durstate[i]);
    printf("\tjumps (out):\n");
    int j = 0;
    do {
      printf("\t\tstate = %d, prob = %f\n",
        seg -> djump_out[i][j], seg -> pjump_out[i][j]);
      j ++;
    } while(seg -> djump_out[i][j - 1] != 1);

    j = 0;
    printf("\tjumps (in):\n");
    do {
      printf("\t\tstate = %d, prob = %f\n",
        seg -> djump_in[i][j], seg -> pjump_in[i][j]);
      j ++;
    } while(seg -> djump_in[i][j - 1] != -1);
  }
}
