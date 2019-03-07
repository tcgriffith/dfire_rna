
GCC = g++
SRCDIR = src
BUILDDIR = .build
BINDIR = bin
INCDIR = include
TARGET = $(BINDIR)/DFIRE_RNA

# ROOT_BOOST = /project/aspen/TC/tools/boost/boost_1_67_0


# SRCEXT = cc
# SOURCES = $(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")

SRCPDB = $(SRCDIR)/PDB.cc \
         $(SRCDIR)/pdb_utils.cc
OBJ_PDB = $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SRCPDB:.cc=.o))

SRCDFIRE = $(SRCDIR)/dfire_calculator.cc \
           $(SRCDIR)/dfire_dihedral.cc \
           $(SRCDIR)/dfire_dipolar.cc \
           $(SRCDIR)/dfire_PDB.cc 
OBJ_DFIRE = $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SRCDFIRE:.cc=.o))

CFLAGS = -std=c++11 -fopenmp# -g -pg
LDFLAGS = -L/home/tc/tools/boost/lib -Llib -fopenmp

LIB = -ldl 

INC = -I include -I $(ROOT_BOOST)


$(TARGET): $(OBJ_PDB) $(OBJ_DFIRE) $(SRCDIR)/main.cc $(INCDIR)/main.h
	@mkdir -p $(BINDIR)
	@echo " Linking..." 
	$(GCC) $^ -o $(TARGET) $(LIB) $(LDFLAGS) $(INC) $(CFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cc $(INCDIR)/%.h
	@mkdir -p $(BUILDDIR)
	$(GCC) $(CFLAGS) $(INC) -c -o $@ $<

.PHONY: clean
clean:
	@echo " Cleaning..."
	rm -rf $(BUILDDIR) $(TARGET) bin/tester

# Tests

# tester: $(OBJ_PDB) $(OBJ_DFIRE) src/tester.cc
# 	$(GCC) $^ $(CFLAGS) $(INC) -static $(LIB) $(LDFLAGS)  -o bin/tester

train: $(OBJ_PDB) $(OBJ_DFIRE) src/train.cc
	$(GCC) $^ $(CFLAGS) $(INC) -static $(LIB) $(LDFLAGS)  -o bin/train

getseq: $(OBJ_PDB) src/getseq.cc
	$(GCC) $^ $(CFLAGS)  $(INC)  $(LIB) $(LDFLAGS)  -o bin/getseq
