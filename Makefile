#CROSS=i686-w64-mingw32-
CC = $(CROSS)gcc
LINK = $(CROSS)gcc
AR = $(CROSS)ar
WITH_SERIALIZATION = true
#CFLAGS = -O1 -g -std=c99 -Wall -fPIC -lm
CFLAGS = -DFP_TYPE=float -Ofast -g -std=c99 -Wall -fPIC -lm $(CFLAGSEXT)
ARFLAGS = -rv
OUT_DIR = ./build
OBJS = $(OUT_DIR)/common.o $(OUT_DIR)/math-funcs.o $(OUT_DIR)/mempool.o $(OUT_DIR)/model.o \
  $(OUT_DIR)/data.o $(OUT_DIR)/generate.o $(OUT_DIR)/inference.o \
  $(OUT_DIR)/estimate.o
ifeq ($(WITH_SERIALIZATION), true)
  OBJS += $(OUT_DIR)/cmp.o $(OUT_DIR)/serial.o
endif
LIBS = -lm

default: $(OUT_DIR)/liblrhsmm.a

test-mempool: test/test-mempool
	test/test-mempool
test: test/test-random-model
	test/test-random-model

$(OUT_DIR)/liblrhsmm.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/liblrhsmm.a $(OBJS)
	@echo Done.

test/test-random-model: test/test-random-model.c $(OUT_DIR)/liblrhsmm.a
	$(LINK) test/test-random-model.c $(OUT_DIR)/liblrhsmm.a $(CFLAGS) \
	  -o test/test-random-model -Wno-unused-function

test/test-mempool: test/test-mempool.c $(OUT_DIR)/liblrhsmm.a
	$(LINK) test/test-mempool.c $(OUT_DIR)/liblrhsmm.a $(CFLAGS) -o test/test-mempool -Wno-unused-function
	test/test-mempool

$(OUT_DIR)/common.o: common.c common.h
$(OUT_DIR)/mempool.o: mempool.c mempool.h
$(OUT_DIR)/model.o: model.c model.h math-funcs.h common.h
$(OUT_DIR)/math-funcs.o: math-funcs.c math-funcs.h
$(OUT_DIR)/data.o: data.c data.h common.h
$(OUT_DIR)/generate.o: generate.c generate.h math-funcs.h common.h
$(OUT_DIR)/inference.o: inference.c inference.h inference-helper.h inference-forward.h \
  inference-forward-geometric.h model.h data.h mempool.h math-funcs.h common.h
$(OUT_DIR)/estimate.o: estimate.c estimate.h inference-helper.h model.h data.h common.h
$(OUT_DIR)/serial.o: serial.c serial.h data.h model.h common.h
$(OUT_DIR)/cmp.o:
	$(CC) $(CFLAGS) -o $(OUT_DIR)/cmp.o -c external/cmp/cmp.c

$(OUT_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $(OUT_DIR)/$*.o -c $*.c

clean:
	@echo 'Removing all temporary binaries...'
	@rm -f $(OUT_DIR)/*.a $(OUT_DIR)/*.o
	@echo Done.

