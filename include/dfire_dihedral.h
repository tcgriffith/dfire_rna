#ifndef _DFIRE_CALC_1
#define _DFIRE_CALC_1
// dfire dihedral term

// #include "PDB.h"
#include "dfire_PDB.h"
#include "pdb_utils.h"
#include "misc.h"
#include "xvec.h"


// #include <unordered_map>
// #include "structure.h"
// #include "xvec.h"
// #include <boost/multi_array.hpp>
// #include <boost/filesystem.hpp>

// ================================
// dfire dihedral term
//
// ================================

class dfire_dihedral{
protected:
    std::string dfire_fpath;
    // void init_dfire_dih();
    // void init_dfire_dih(std::string);

    void init_dfire_etable(std::string);

    // void set_dih_reftype(int reftype){dih_ref_type = reftype;}
    // int get_dih_reftype(){return dih_ref_type;}

    static constexpr int dbin=80, dtype=7;
    static constexpr double dbin_width=4.5;
    int dih_ref_type = 0;

    int count_dih[dtype][dbin]={};
    double dfire_dih[dtype][dbin]={};

    std::unordered_map<std::string,int> dihtype;
    std::unordered_map<int,std::string> dihtype2;

    double dihedral(atom &aa, atom &ab, atom &ac, atom &ad);

// dihedral term

    // reference type: 0 -- average, as in 3drnascore
    //                 1 -- uniform distribution
    // used in count_to_dfiredih



    // void count_dihs(chain & achain); // deprecated
    // void count_dihs(structure & astr);
    void count_dihs(rna & astr);
    void count_to_dfiredih();

    // uniform distribution
    // void count_to_dfiredih_unibg();

    int get_dihtype_int(std::string);
    // void init_dfire_dih(std::string);


    // double calc_dfire_dih(chain & achain); // deprecated

    // double calc_dfire_dih(structure & astr);
    double calc_dfire_dih(rna & astr);



public:

    dfire_dihedral(){};
    dfire_dihedral(std::string energy_dir){
        init_dfire_etable(energy_dir);
    };

    ~dfire_dihedral(){};
    // ==============================================
    // training
    // ==============================================

    void train_dfire_dih(string pdbl_fpath, int length_limit = 500);
    void print_edfire_dih(FILE * out);
    void print_count_dih(FILE * out);

    // ==============================================
    // energy calculation
    // input: structure, may have multiple chains
    // output: dfire dihedral term score
    // ==============================================


    // double get_dfire(structure & astr){return calc_dfire_dih(astr);};
    double get_dfire(rna & astr){return calc_dfire_dih(astr);};
    void get_dfire_per_res(structure & astr, FILE *out);
};



#endif
