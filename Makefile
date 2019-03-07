
GCC = g++
SRCDIR = src
BUILDDIR = .build
BINDIR = bin
INCDIR = include
TARGET = $(BINDIR)/DFIRE_RNA

SRCPDB = $(SRCDIR)/PDB.cc \
         $(SRCDIR)/pdb_utils.cc
OBJ_PDB = $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SRCPDB:.cc=.o))

SRCDFIRE = $(SRCDIR)/dfire_calculator.cc \
           $(SRCDIR)/dfire_PDB.cc 
OBJ_DFIRE = $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SRCDFIRE:.cc=.o))

CFLAGS = -std=c++11 -fopenmp# -g -pg
LDFLAGS = -Llib -fopenmp

LIB = -ldl 

INC = -I include


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
	rm -rf $(BUILDDIR) $(BINDIR)

# For training
train: $(OBJ_PDB) $(OBJ_DFIRE) src/train.cc
	@mkdir -p $(BINDIR)
	$(GCC) $^ $(CFLAGS) $(INC) -static $(LIB) $(LDFLAGS)  -o bin/train

# get fasta sequences from PDB
getseq: $(OBJ_PDB) src/getseq.cc
	@mkdir -p $(BINDIR)
	$(GCC) $^ $(CFLAGS)  $(INC)  $(LIB) $(LDFLAGS)  -o bin/getseq
