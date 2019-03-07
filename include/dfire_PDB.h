#ifndef _DFIRE_PDB
#define _DFIRE_PDB

#include "PDB.h"

// class atom;
class rnabase;
class rnachain;
class rna;



class rnabase: public residue{
protected:
    bool _has_dih=false;
    bool _miss_atom=false;
    bool _has_pair=false;

    double e_atom = 0.0;
    double e_dih = 0.0;
    double e_dipo = 0.0;

        // double alpha, beta, gamma, delta, epsilon, zeta, chi;
    double dihedrals[7] = {0,0,0,0,0,0,0};
public:
    // std::vector<atom> atoms;
    rnabase(){};
    ~rnabase(){};

    bool has_dih(){return _has_dih;};
    
    double * get_dihedrals(){
        return dihedrals;
    }

    inline std::string dihedrals_to_string(){
        // for (int i = 0; i < 7; ++i)

        char tmpstring[500];

            sprintf(tmpstring,"%s %d %s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", resname.c_str(), _has_dih, resid.c_str(),
            dihedrals[0],dihedrals[1],dihedrals[2],dihedrals[3],dihedrals[4],dihedrals[5],dihedrals[6]);
        return string(tmpstring);
    }

    void calc_torsion(rnabase *r_up, rnabase * r_down);

};

class rnachain: public chain{
public:
    std::vector<rnabase> residues;
    rnachain(){};
    ~rnachain(){};

    bool init_dihedrals();

};

class rna: public structure{
public:
    std::vector<rnachain> chains;
    int n_atompairs = 0;
    rna(){};
    rna(std::string pdbfile){
        readpdb(pdbfile);
    }
    ~rna(){};


    bool readpdb(std::string);

    bool rna_init_polar();

    bool rna_init_pairs();

    bool rna_init_dihedrals() { for (auto & achain : chains) achain.init_dihedrals(); return true;};


    // bool rna_init_dihedrals();
    // bool rna_init_dihedrals() { for (auto & achain : chains) achain.init_dihedrals(); return true;};
};


#endif