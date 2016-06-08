%{
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
%}

printf("#include \"math-funcs.h\"\n\n");
printf("FP_TYPE lrh_lse_table[lrh_lse_table_size] = {\n  ");

tbsize = 7000;
tbstep = 0.001;

x = -(0:tbsize - 1) * tbstep;
y = log(1.0 + exp(x));

for n = 1:tbsize
  if(mod(n, 5) == 0)
    if(n == tbsize)
      printf("%2.7e\n};\n", y(n));
    else
      printf("%2.7e,\n  ", y(n));
    end
  else
    printf("%2.7e, ", y(n));
  end
end

printf("\n");

tbsize = 7000;
tbstep = 0.002;

x = -(0:tbsize - 1) * tbstep;
y = exp(x);

printf("FP_TYPE lrh_exp_table[lrh_exp_table_size] = {\n  ");
for n = 1:tbsize
  if(mod(n, 5) == 0)
    if(n == tbsize)
      printf("%2.7e\n};\n", y(n));
    else
      printf("%2.7e,\n  ", y(n));
    end
  else
    printf("%2.7e, ", y(n));
  end
end

printf("\n");

