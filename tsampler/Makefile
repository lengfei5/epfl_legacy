PREFIX = /mnt/common/epfl
SAMLIBDIR = $(PREFIX)/lib
SAMINCLUDE = $(PREFIX)/include
CXXFLAGS = -pg -O3 -Wall -I$(SAMINCLUDE) 
LIBDIRS = -L$(SAMLIBDIR)
LDLIBS = $(LIBDIRS) -lbam -lz -lgsl -lgslcblas

OBJ = gene.o tokens.o main.o
BIN = rna

all: $(BIN)

rna: $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDLIBS)

install: all
	install $(BIN) $(PREFIX)/bin

clean:
	$(RM) $(OBJ) $(BIN)

cleaner: clean
	$(RM) core *~
