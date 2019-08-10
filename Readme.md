

# DFIRE-RNA

Knowledge-based statistical energy based on DFIRE reference state for RNA structure evaluation.

# Install

Use cmake

```
# cmake stuff

mkdir build
cd build
cmake ..
make

# set DFIRE_RNA_HOME in .bashrc (source it or reopen the terminal)
./install.sh

```

## For Linux

Use Makefile (use gcc compiler):

```
make
# You can manually set environment variable DFIRE_RNA_HOME to the source directory of DFIRE-RNA so the energy files can be detected, or run 
```

# Usage

```
# Help is displayed if run dfire_rna without parameters

#######################################################
# Calculate dfire_rna score for a pdb or a list of pdbs
#######################################################
Usage: ./bin/dfire_rna pdb 
   or: ./bin/dfire_rna [ options ] 
Options:
   pdb [ pdb2 pdb3 ...], input RNA structures in pdb format
   -d directory,         OPTIONAL, override default directory of energyfiles
                         default: dfire_rna/data/energyfiles
   -norm                 normalize DFIRE score by RNA length (experimental)
   -l pdblist,           A list of absolute paths to pdb files (plain text) UTF-8 encoding

```


# Example of usage

```
# Run in the example dir

../bin/DFIRE_RNA *.pdb

# output:
#
# 1a9nR.pdb -12480.188898  
# 3b58ABC.pdb -53013.816472  
# 5e3hBC.pdb -19451.089802  
# S_000001_000.pdb -52615.344214  

```




# Citation

T. Zhang, G. Hu, Y. Yang, J. Wang, and Y. Zhou, “All-atom knowledge-based potential for RNA structure discrimination based on the distance-scaled finite ideal-gas reference state.”, J. Computational Biology, in press (2019).

More works related to DFIRE reference state can be found at:

http://sparks-lab.org/index.php/Main/Publications