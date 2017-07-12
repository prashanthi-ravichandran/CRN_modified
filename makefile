# makefile for CardioWave
#
#	created by cw_make.pl
#
#
# compiler infomation
CC =  gcc
CFLAGS =  -O5 -I/usr/include/malloc
LD =  gcc
LFLAGS =  -O5 -I/usr/include/malloc
#
# targets
HDR =	CardioWave.h VectorOps.h
OBJ =	MemCRN.o Domain123.o TstepEMrk23_tol.o \
	StimTrain.o OutputDump.o ParallelNone.o \
	DebugNone.o CWaveKernel.o VectorOps.o \
	StencilOps1D.o 
LIB =	 -lm
#
# the main dependency list
Run: $(HDR) Run.c $(OBJ)
	$(LD) $(LFLAGS) -o $@ Run.c $(OBJ) $(LIB)
	@echo ____DONE____
#
# create the Run.c file with cw_run.pl
Run.c:
	perl $(CWAVEHOME)/scripts/cw_run.pl $(OBJ)
CWaveKernel.o: CWaveKernel.c
	$(CC) $(CFLAGS) -o $@ -DCWAVEARCH_Darwin -c $<
#
# create the sample input file
sample.in:
	perl $(CWAVEHOME)/scripts/cw_sample.pl $(OBJ)
#
# default .c to .o compilation
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
#
# clean up
clean:
	@perl $(CWAVEHOME)/scripts/cw_clean.pl $(OBJ)
	@echo _ALL_CLEAN_
#
# compare files back to CWAVE repository
compare:
	@perl $(CWAVEHOME)/scripts/cw_cmp.pl $(OBJ)
#
# get files from the archives
getfiles:
	@perl $(CWAVEHOME)/scripts/cw_get.pl $(OBJ)
	@echo _FILES_RETRIEVED_

