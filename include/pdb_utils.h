#ifndef _PDB_UTILS
#define _PDB_UTILS
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>


// using namespace std;

// Necessary parameters for RNA structures

namespace pdb_utils{

    extern double hbond_cutoff; // hydrogen bond cutoff

// ==============================================
// smcra/RNA
// special RNA NTs
// ==============================================

   extern std::unordered_map<std::string,std::string> rna_specials;


// ==============================================
// smcra/RNA
// Deal with old/new PDB conventions
// ==============================================

   extern std::unordered_map<std::string,std::string> old_2_new;

// ==============================================
// neighbor atom maps in RNA
// ==============================================

   extern std::unordered_map<std::string, std::vector<std::string>> rna_atom_neighbors;


// ==============================================
// Protein Triple code to single char code
// ==============================================


  extern std::unordered_map<std::string, std::string> prot_tri2sin;

//  ==============================================
//  RNA polar atom type encoding
//  RNA {atom name,int}  <-> {int,atom name}
//             polartype <-> typepolar
//  ==============================================

  extern std::unordered_map<std::string, int> polartype;

  extern std::unordered_map<int, std::string> typepolar;


//  ==============================================
//  RNA atom type encoding
//  RNA {atom name,int} <-> {int,atom name}
//             atomtype <-> typeatom
//  ==============================================
  extern std::unordered_map<std::string,int> atomtype;

  extern std::unordered_map<int, std::string> typeatom;

//  ==============================================
//  functions
//  ==============================================

    
    // 0 Protein
    // 1 RNA
    // 2 other
    int guess_ctype(const std::string& resname);

    bool is_rna_polar_atom(const std::string& atype);

    bool is_rna_sideatom(const std::string& atype_s);

    int get_polartype_int(std::string atype_s);

    int get_atomtype_int(std::string atype_s);
    // inline int get_atomtype_int(std::string atype_s){
    //     if (atomtype.find(atype_s) != atomtype.end()){

    //         return atomtype[atype_s];
    //     }
    //     else{
    //         std::cerr << atype_s << " not found" << std::endl;
    //         return -1;
    //     }

    // }

    

// pdb format conversion, eg C4* -> C4', see old_2_new
    std::string format_rna(const std::string oldstr);

    std::string protein_tri2sin(const std::string &str_tri);

//  calculate pdb
    double calc_RMSD();

}

#endif
