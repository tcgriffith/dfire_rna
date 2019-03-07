#ifndef _DFIRE_CALC_2
#define _DFIRE_CALC_2

#include <unordered_map>

// #include "PDB.cc"
#include "dfire_PDB.h"
#include "pdb_utils.h"
#include "misc.h"
#include "xvec.h"

// #include "chain.h"
// #include "dfire_calculator.h"
// #include "structure.h"
// #include "xvec.h"
#include <boost/multi_array.hpp>
// #include <boost/filesystem.hpp>




// ================================
// dfire dipolar atom term
//
// ================================

class dfire_dipolar{

public:
// ==============================================
// dfire dipole
// ==============================================

    // smcra_param * param_p;
    static constexpr double ALPHA=1.61;
    static constexpr double PI=3.1415926, radian=180.0/3.1415926;
    static constexpr int dipo_bin=6, dipo_type=43;
    static constexpr double dipo_width = 1.0/3;
    static constexpr int polartypen = 39, matype = 85, mbin = 30;



// dipolar bin divided by cos(theta)

// c++ complains about stack overflow for the following multi arrays,
// so boost::multi_array is used here

    // double dfire_dipo_pn[matype][matype][dipo_bin][mbin]={};

    // int count_dipo_pp[matype][matype][dipo_bin][mbin]={};
    // double dfire_dipo_pp[matype][matype][dipo_bin][mbin]={};

    typedef boost::multi_array<int, 4> int_array;
    typedef boost::multi_array<double, 4> double_array;
    typedef int_array::index index_i;
    typedef double_array::index index_d;
// polar - nonpolar
    int_array count_dipo_pn{boost::extents[polartypen][matype][dipo_bin][mbin]};
    double_array dfire_dipo_pn{boost::extents[polartypen][matype][dipo_bin][mbin]};
// polar -polar
    int_array count_dipo_pp{boost::extents[polartypen][polartypen][dipo_bin][mbin]};
    double_array dfire_dipo_pp{boost::extents[polartypen][polartypen][dipo_bin][mbin]};

    // typedef std::vector<std::vector<std::vector<std::vector<int >>>> int_array;
    // typedef std::vector<std::vector<std::vector<std::vector<double >>>> double_array;



    void init_dipo_arrays();
    void init_dfire_etable_dipo(std::string);


// dipo
    void train_dfire_dipo(string pdbl_fpath, int length_limit = 1000);

    // deprecated
    // double calc_dfire_dipo(chain & achain);

    std::pair<double, double> calc_dfire_dipo(rna & astr);
    void init_dfire_dipo();
    void init_dfire_dipo(std::string);

    void print_count_dipo(string &,string &);
    void print_edfire_dipo(string &, string &);

    void print_count_dipo(FILE * out_pn, FILE * out_pp);
    void print_edfire_dipo(FILE * out_pn, FILE * out_pp);

    // deprecated
    // void count_dipo(chain &);

    // void count_dipo(structure & astr);
    void count_dipo(rna &astr);
    void count_to_dfiredipo();

    int ang2bin(double angle);
    int cosang2bin(double);



// ==============================================
// dfire dipolar
// use get_dfire(structure & astr) to get energy
// ==============================================

    dfire_dipolar(std::string energy_dir){
        // param_p = new smcra_param();
        init_dfire_etable_dipo(energy_dir);
    }


    dfire_dipolar(){};
    ~dfire_dipolar(){};

    // double get_dfire(structure & astr){return calc_dfire_dipo(astr);};
    // std::pair<double, double> get_dfire(structure & astr){return calc_dfire_dipo(astr);};
   std::pair<double, double> get_dfire(rna & astr){return calc_dfire_dipo(astr);};
    void get_dfire_per_res(structure & astr, FILE *out);

    void get_dfire_profile(structure & astr, FILE *out);
};

#endif
