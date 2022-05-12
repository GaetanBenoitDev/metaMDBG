#CC          = gcc
EXTRA_FLAGS = -Wno-unused-function -Wno-misleading-indentation -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES
CFLAGS      = -Wall -O3 $(EXTRA_FLAGS)

SIMD_FLAG   = -march=native

ifneq ($(armv7),) # for ARMv7
	SIMD_FLAG   =  -march=armv7-a -mfpu=neon -D__AVX2__
else
ifneq ($(armv8),) # for ARMv8
ifneq ($(aarch64),) # for Aarch64 
	SIMD_FLAG   =  -march=armv8-a+simd -D__AVX2__
else # for Aarch32
	SIMD_FLAG   =  -march=armv8-a+simd -mfpu=auto -D__AVX2__
endif
endif
endif

# for debug
ifneq ($(debug),)
	DFLAGS   =   -D __DEBUG__
endif
# for gdb
ifneq ($(gdb),)
	CFLAGS   = -Wall -g ${DFLAGS} $(EXTRA_FLAGS)
else
	CFLAGS   = -Wall -O3 ${DFLAGS} $(EXTRA_FLAGS)
endif

# for gprof
ifneq ($(pg),)
	PG_FLAG  =   -pg
	CFLAGS  +=   -pg
endif

LIB     = -lm -lz -lpthread
ifneq ($(PREFIX),)
	OUT_PRE_DIR = $(PREFIX)
else
	OUT_PRE_DIR = .
endif

BIN_DIR = $(OUT_PRE_DIR)/bin
LIB_DIR = $(OUT_PRE_DIR)/lib
INC_DIR = ./include
SRC_DIR = ./src

SOURCE = $(SRC_DIR)/abpoa_align.c $(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa_graph.c $(SRC_DIR)/abpoa_plot.c $(SRC_DIR)/abpoa_seed.c $(SRC_DIR)/abpoa_seq.c $(SRC_DIR)/abpoa_output.c $(SRC_DIR)/kalloc.c $(SRC_DIR)/kstring.c  $(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/simd_check.c $(SRC_DIR)/utils.c
HEADER = $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa.h $(INC_DIR)/abpoa.h $(SRC_DIR)/abpoa_seed.h $(SRC_DIR)/abpoa_seq.h $(SRC_DIR)/abpoa_output.h $(SRC_DIR)/kalloc.h $(SRC_DIR)/kdq.h $(SRC_DIR)/khash.h $(SRC_DIR)/kseq.h $(SRC_DIR)/ksort.h $(SRC_DIR)/kstring.h $(SRC_DIR)/kvec.h $(SRC_DIR)/simd_instruction.h $(INC_DIR)/simd_instruction.h $(SRC_DIR)/simd_abpoa_align.h $(SRC_DIR)/utils.h
OBJS   = $(SRC_DIR)/abpoa_align.o $(SRC_DIR)/abpoa_graph.o $(SRC_DIR)/abpoa_plot.o $(SRC_DIR)/abpoa_seed.o $(SRC_DIR)/abpoa_seq.o $(SRC_DIR)/abpoa_output.o $(SRC_DIR)/kalloc.o $(SRC_DIR)/kstring.o $(SRC_DIR)/simd_abpoa_align.o $(SRC_DIR)/simd_check.o $(SRC_DIR)/utils.o

# SIMD label
SIMD_CHECK_D = -D __CHECK_SIMD_MAIN__

FLAG_SSE2     = -msse2
FLAG_SSE41    = -msse4.1
FLAG_AVX2     = -mavx2
# FLAG_AVX512F  = -mavx512f
# FLAG_AVX512BW = -mavx512bw

ifneq ($(sse2),)
	SIMD_FLAG=$(FLAG_SSE2)
	py_SIMD_FLAG = SSE2=1
else ifneq ($(sse41),)
	SIMD_FLAG=$(FLAG_SSE41)
	py_SIMD_FLAG = SSE41=1
else ifneq ($(avx2),)
	SIMD_FLAG=$(FLAG_AVX2)
	py_SIMD_FLAG = AVX2=1
#else ifneq ($(avx512f),)
#	SIMD_FLAG=$(FLAG_AVX512F)
#	py_SIMD_FLAG = AVX512f=1
#else ifneq ($(avx512bw),)
#	SIMD_FLAG=$(FLAG_AVX512BW)
#	py_SIMD_FLAG = AVX512BW=1
endif

.c.o:
		$(CC) -c $(CFLAGS) $< -I$(INC_DIR) -o $@

BIN      = $(BIN_DIR)/abpoa
ifneq ($(gdb),)
	BIN  = $(BIN_DIR)/gdb_abpoa
endif
ABPOALIB = $(LIB_DIR)/libabpoa.a
# TODO add example
EXAMPLE  = example


all:       $(BIN) 
abpoa:     $(BIN)
libabpoa:  $(ABPOALIB)
example:   $(EXAMPLE)

$(BIN):$(SRC_DIR)/abpoa.o $(ABPOALIB)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(CFLAGS) $< -I$(INC_DIR) -L$(LIB_DIR) -labpoa $(LIB) -o $@ $(PG_FLAG)

$(EXAMPLE):example.c $(ABPOALIB)
	$(CC) $(CFLAGS) $< -o $@ -I$(INC_DIR) -L$(LIB_DIR) -labpoa $(LIB)

$(ABPOALIB):$(OBJS)
	if [ ! -d $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi
	$(AR) -csr $@ $(OBJS)

$(SRC_DIR)/abpoa.o:$(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h \
                   $(SRC_DIR)/abpoa_seq.h $(SRC_DIR)/utils.h $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) -I$(INC_DIR) $< -o $@

$(SRC_DIR)/simd_check.o:$(SRC_DIR)/simd_check.c $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) -I$(INC_DIR) $< -o $@

$(SRC_DIR)/simd_abpoa_align.o:$(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/simd_instruction.h $(SRC_DIR)/utils.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) -I$(INC_DIR) $< -o $@

install_py: python/cabpoa.pxd python/pyabpoa.pyx python/README.md
	${py_SIMD_FLAG} python setup.py install
	
sdist: install_py
	${py_SIMD_FLAG} python setup.py sdist #bdist_wheel

publish_pypi: clean_py sdist
	twine upload dist/*

clean:
	rm -f $(SRC_DIR)/*.[oa] $(LIB_DIR)/*.[oa] $(BIN)
clean_py:
	rm -rf build/ dist/ pyabpoa.egg-info/ python/pyabpoa.c
