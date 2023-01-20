IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -O3 -std=c++17 -Wall -D SIZE=$(SIZE) -lpthread

ODIR=../obj

_DEPS = skinny.h trails.h keyDependent.h printDetails.h experiment.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main$(SIZE).o trails$(SIZE).o skinny$(SIZE).o keyDependent$(SIZE).o  printDetails$(SIZE).o experiment$(SIZE).o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BDIR=../bin



# ${BDIR}/test: test.cpp $(OBJ)
# 	$(CC) -o $@ $^ $(CFLAGS)

# ${BDIR}/experiment$(SIZE): experiment.cpp $(ODIR)/skinny$(SIZE).o $(ODIR)/trails$(SIZE).o $(ODIR)/keyDependent$(SIZE).o $(ODIR)/printDetails$(SIZE).o
# 	$(CC) -o $@ $^ $(CFLAGS)

${BDIR}/skinnyKeyDependent$(SIZE): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

${ODIR}/main$(SIZE).o: main.cpp  ${ODIR}/printDetails$(SIZE).o ${ODIR}/skinny$(SIZE).o ${ODIR}/keyDependent$(SIZE).o ${ODIR}/experiment$(SIZE).o
	$(CC) -c -o $@ $< $(CFLAGS)

${ODIR}/experiment$(SIZE).o: experiment.cpp ${IDIR}/experiment.h ${ODIR}/keyDependent$(SIZE).o ${ODIR}/trails$(SIZE).o ${ODIR}/skinny$(SIZE).o ${ODIR}/printDetails$(SIZE).o
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/keyDependent$(SIZE).o: keyDependent.cpp ${IDIR}/keyDependent.h ${ODIR}/trails$(SIZE).o  ${ODIR}/skinny$(SIZE).o
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/trails$(SIZE).o: trails.cpp ${IDIR}/trails.h ${ODIR}/skinny$(SIZE).o
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/skinny$(SIZE).o: skinny.cpp ${IDIR}/skinny.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/printDetails$(SIZE).o: printDetails.cpp ${IDIR}/printDetails.h ${ODIR}/skinny$(SIZE).o ${ODIR}/keyDependent$(SIZE).o
	$(CC) -c -o $@ $< $(CFLAGS)




.PHONY: clean

clean:
	rm -f $(BDIR)/*  $(ODIR)/*.o *~ core $(IDIR)/*~ 