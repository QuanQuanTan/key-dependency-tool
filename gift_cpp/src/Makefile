IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -O3 -std=c++17 -D SIZE=$(SIZE) -lcryptominisat5

ODIR=../obj

_DEPS = gift.h trails.h keyDependent.h gaussianElimination.h SAT.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = gift$(SIZE).o trails$(SIZE).o gaussianElimination${SIZE}.o SAT${SIZE}.o keyDependent$(SIZE).o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BDIR=../bin

# ${BDIR}/test: test.cpp $(OBJ)
# 	$(CC) -o $@ $^ $(CFLAGS)

${BDIR}/giftKeyDependent$(SIZE): main.cpp $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

$(ODIR)/keyDependent$(SIZE).o: keyDependent.cpp ${IDIR}/keyDependent.h ${IDIR}/trails.h ${IDIR}/gift.h ${IDIR}/gaussianElimination.h ${IDIR}/SAT.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/gaussianElimination${SIZE}.o: gaussianElimination.cpp ${IDIR}/gaussianElimination.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/SAT$(SIZE).o: SAT.cpp ${IDIR}/SAT.h ${IDIR}/gift.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/trails$(SIZE).o: trails.cpp ${IDIR}/trails.h ${IDIR}/gift.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/gift$(SIZE).o: gift.cpp ${IDIR}/gift.h
	$(CC) -c -o $@ $< $(CFLAGS)




.PHONY: clean

clean:
	rm -f $(BDIR)/*  $(ODIR)/*.o *~ core $(IDIR)/*~ 