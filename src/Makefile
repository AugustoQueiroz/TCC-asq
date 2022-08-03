IDIR = ../include
CC=gcc
CFLAGS=-I$(IDIR) -Wall -fopenmp

BUILDDIR=build

LIBS=-lm

_DEPS = CountMin.h Hashing.h KMerProcessing.h SequencingSimulator.h mathutils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_MAINOBJ = main.o CountMin.o Hashing.o KMerProcessing.o mathutils.o
MAINOBJ = $(patsubst %,$(BUILDDIR)/%,$(_MAINOBJ))

_SIMULATOROBJ = ReadsGenerator.o SequencingSimulator.o KMerProcessing.o
SIMULATOROBJ = $(patsubst %,$(BUILDDIR)/%,$(_SIMULATOROBJ))

_TRAVERSALOBJ = CountMin.o Hashing.o KMerProcessing.o mathutils.o
TRAVERSALOBJ = $(patsubst %,$(BUILDDIR)/%,$(_TRAVERSALOBJ))

$(BUILDDIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(MAINOBJ)
	$(CC) -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

simulator: $(SIMULATOROBJ)
	$(CC) -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

traversal: traversal.cpp $(TRAVERSALOBJ)
	g++ -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

counting: counting.cpp $(TRAVERSALOBJ)
	g++ -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

construct-hashtable: ConstructDBHTFromDBCM.cpp $(TRAVERSALOBJ) HashTable.o
	g++ -o $(BUILDDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm $(BUILDDIR)/*.o
