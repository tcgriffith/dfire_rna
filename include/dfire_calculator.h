#ifndef _DFIRE_CALC
#define _DFIRE_CALC

#include <unordered_map>
// #include "structure.h"
// #include "PDB.h"
#include "dfire_PDB.h"
#include "pdb_utils.h" 
#include "misc.h"
#include "xvec.h"

#include <fstream>
// #include <boost/filesystem.hpp>
// #include "SMCRA_param.h" // already included in structure.h
// #include "boost/multi_array.hpp" // moved to dfire_dipolar.h

// ========================
// vanilla dfire atom-atom distance term for rna
//
// ========================

// #define DEBUG

class dfire_calculator{
protected:
// ==============================================
// common

    // std::string dfire_fpath;

    // deprecated, moved to SMCRA_param.h, use param_p
    // std::unordered_map<std::string,int> atomtype;
    // std::unordered_map<int,std::string> atomtype2;
    // void init_rna_atomtype();
    // int get_atomtype_int(std::string);
    // smcra_param * param_p;

// ==============================================
// dfire vanilla
// params and etable
// modified: 181115 
// make alpha, bin_size, max_dist adjustable for better tweaking.
// ==============================================
     
    static constexpr int mbin=100, matype = 85;


    double ALPHA = 1.61;
    double bin_size = 0.5;
    double max_dist = 15.0;
    int max_bin = 30;




    double dfire[matype][matype][mbin]={};
    int count[matype][matype][mbin]={};
    double dfire_c[matype][matype][mbin]={};
// ==============================================
// dfire vanilla
// init
// ==============================================

    // void init_dfire();
    // void init_dfire(std::string);
    // void init_dfire_etable();
    


// ==============================================
// dfire vanilla
// training
// ==============================================


    
    void count_to_dfire();
    // void count_pairs(structure &);
    void count_pairs(rna &);

    bool is_far_away(atom & atom_a, atom & atom_b, double maxdist);

    // 
    // void count_to_3drna();
    // double calc_3drna();

    // double calc_dfire(chain & achain); // deprecated
    // double calc_dfire(structure & astr); // deprecated, use dfire_PDB::rna instead;

    double calc_dfire(rna &astr);





public:

// ==============================================
// dfire vanilla
// use get_dfire(structure & astr) to get dfire energy
// ==============================================

    // used for training
    // dfire_calculator(){param_p = new smcra_param();};
    dfire_calculator(){};


    dfire_calculator(const std::string &energy_dir){
        // param_p = new smcra_param();
        init_dfire_etable(energy_dir);
    };
    void init_dfire_etable(const std::string &energy_dir);




    ~dfire_calculator(){};
    void train_dfire(string pdbl_fpath, double alpha, double bin_size, double max_dist);

// alpha, bin_size, max_dist
// modified: 181115 
    void set_abm(double param_alpha, double param_bin_size, double param_max_dist){
        this->ALPHA = param_alpha;
        this->bin_size = param_bin_size;

        double residual = fmod(param_max_dist, param_bin_size);
        int mm = (int) (param_max_dist / param_bin_size);

        // adjust the last bin to make it a complete bin
        if (residual > 1e-6) {
            cerr << "# max_dist adjusted to upperbound" << endl;
            mm ++;
        }

        this->max_bin = mm;
        this->max_dist = mm * bin_size;
    }

// translate distance to bin, based on bin_size
// modified: 181115 
    inline int distance2bin(double dist, double bin_size){
        int bin = (int) (dist / bin_size);
        return bin;
    }

    void print_parameters(FILE * out);

    void print_edfire(FILE *out);
    void print_count(FILE *out);


    // double get_dfire(structure & astr){return calc_dfire(astr);}; //deprecated, use dfire_PDB
    double get_dfire(rna & astr){return calc_dfire(astr);};
    void get_dfire_per_res(structure & astr, FILE *out);
    

    // similar to 3drna;
    void train_3drna(string pdbl_fpath);
    void count_to_3drna();
    double  calc_3drna(rna &astr);
    double get_3drna(rna & astr){return calc_3drna(astr);};



};

#endif
