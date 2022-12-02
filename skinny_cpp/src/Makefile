IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -O3 -std=c++17 -D SIZE=$(SIZE)

ODIR=../obj

_DEPS = skinny.h trails.h keyDependent.h printDetails.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = skinny$(SIZE).o trails$(SIZE).o keyDependent$(SIZE).o printDetails.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BDIR=../bin



# ${BDIR}/test: test.cpp $(OBJ)
# 	$(CC) -o $@ $^ $(CFLAGS)

${BDIR}/skinnyKeyDependent${SIZE}: main.cpp $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)



$(ODIR)/keyDependent$(SIZE).o: keyDependent.cpp ${IDIR}/keyDependent.h ${IDIR}/trails.h ${IDIR}/skinny.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/trails$(SIZE).o: trails.cpp ${IDIR}/trails.h ${IDIR}/skinny.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/skinny$(SIZE).o: skinny.cpp ${IDIR}/skinny.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/printDetails.o: printDetails.cpp ${IDIR}/printDetails.h ${IDIR}/skinny.h ${IDIR}/keyDependent.h 
	$(CC) -c -o $@ $< $(CFLAGS)




.PHONY: clean

clean:
	rm -f $(BDIR)/*  $(ODIR)/*.o *~ core $(IDIR)/*~ 