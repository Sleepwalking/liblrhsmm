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
#include "mempool.h"

static lrh_mempool_page* lrh_create_mempool_page(size_t capacity) {
  lrh_mempool_page* ret = malloc(sizeof(lrh_mempool_page));
  ret -> data = malloc(capacity);
  ret -> top = 0;
  ret -> capacity = capacity;
  ret -> next = NULL;
  return ret;
}

lrh_mempool* lrh_create_mempool(size_t pagesize) {
  lrh_mempool* ret = malloc(sizeof(lrh_mempool));
  ret -> pagesize = pagesize;
  ret -> totalsize = pagesize;
  ret -> npage = 1;
  ret -> first = lrh_create_mempool_page(pagesize);
  ret -> last = ret -> first;
  return ret;
}

void lrh_delete_mempool(lrh_mempool* dst) {
  if(dst == NULL) return;
  lrh_mempool_page* curr = dst -> first;
  do {
    free(curr -> data);
    lrh_mempool_page* next = curr -> next;
    free(curr);
    curr = next;
  } while(curr != NULL);
  free(dst);
}

void* lrh_mempool_malloc(lrh_mempool* src, size_t size) {
  void* ret = NULL;
  if(src -> last -> top + size <= src -> last -> capacity) {
    ret = (char*)src -> last -> data + src -> last -> top;
    src -> last -> top += size;
  } else {
    lrh_mempool_page* newpage = NULL;
    if(size > src -> pagesize)
      newpage = lrh_create_mempool_page(size);
    else
      newpage = lrh_create_mempool_page(src -> pagesize);
    src -> last -> next = newpage;
    src -> last = newpage;
    src -> npage ++;
    src -> totalsize += newpage -> capacity;
    ret = newpage -> data;
    newpage -> top += size;
  }
  return ret;
}

void* lrh_mempool_calloc(lrh_mempool* src, size_t number, size_t size) {
  void* ret = lrh_mempool_malloc(src, size * number);
  memset(ret, 0, size * number);
  return ret;
}

