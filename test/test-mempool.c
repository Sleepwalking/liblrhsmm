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
#include <stdlib.h>
#include "../mempool.h"

int main() {
  lrh_mempool* pool = lrh_create_mempool(1024 * 1024); //1Mb page size
  for(int i = 0; i < 1e7; i ++) { // 10 million mallocs
    int isize = rand() % 50; // allow zero size malloc
    if(i == 1000 || i == 100000)
      isize = 1024 * 1024 * 5; // 2Mb * 2 calloc request
    char* tmp = lrh_mempool_calloc(pool, isize, 2);
    if(isize > 0)
      tmp[0] = 0;
  }
  printf("Total size = %ld bytes.\n", pool -> totalsize);
  printf("Number of pages = %ld.\n", pool -> npage);
  lrh_delete_mempool(pool);
  return 0;
}

