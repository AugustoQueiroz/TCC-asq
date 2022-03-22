IDIR = ../include
CC=gcc
CFLAGS=-I$(IDIR) -Wall

BUILDDIR=build

LIBS=-lm

_DEPS = CountMin.h Hashing.h KMerProcessing.h SequencingSimulator.h mathutils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o CountMin.o Hashing.o KMerProcessing.o SequencingSimulator.o mathutils.o
OBJ = $(patsubst %,$(BUILDDIR)/%,$(_OBJ))

$(BUILDDIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

simulator: $(BUILDDIR)/ReadsGenerator.o $(BUILDDIR)/SequencingSimulator.o
	$(CC) -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm $(BUILDDIR)/*.o