IDIR = ../include
CC=gcc
CFLAGS=-I$(IDIR) -Wall

BUILDDIR=build

LIBS=-lm

_DEPS = CountMin.h Hashing.h KMerProcessing.h SequencingSimulator.h mathutils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_MAINOBJ = main.o CountMin.o Hashing.o KMerProcessing.o mathutils.o
MAINOBJ = $(patsubst %,$(BUILDDIR)/%,$(_MAINOBJ))

_SIMULATOROBJ = ReadsGenerator.o SequencingSimulator.o
SIMULATOROBJ = $(patsubst %,$(BUILDDIR)/%,$(_SIMULATOROBJ))

$(BUILDDIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(MAINOBJ)
	$(CC) -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

simulator: $(SIMULATOROBJ)
	$(CC) -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm $(BUILDDIR)/*.o