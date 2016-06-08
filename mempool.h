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

#ifndef LRHSMM_MEMPOOL_H
#define LRHSMM_MEMPOOL_H

#include <unistd.h>

// linked list based memory pool implementation for reducing time spent on
//   frequent mallocs in the inference of state duration and mixture occupancy
//   probabilities

typedef struct lrh_mempool_page {
  struct lrh_mempool_page* next;
  void* data;
  size_t top;
  size_t capacity;
} lrh_mempool_page;

typedef struct {
  lrh_mempool_page* first;
  lrh_mempool_page* last;
  size_t pagesize;
  size_t totalsize;
  size_t npage;
} lrh_mempool;

lrh_mempool* lrh_create_mempool(size_t pagesize);
void lrh_delete_mempool(lrh_mempool* dst);

void* lrh_mempool_malloc(lrh_mempool* src, size_t size);
void* lrh_mempool_calloc(lrh_mempool* src, size_t number, size_t size);

#endif

