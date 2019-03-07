// #include "structure.h"
// #include <boost/filesystem.hpp>

#include "PDB.h"
#include "dfire_PDB.h"
#include <algorithm>
#include <unordered_set>
#include <cstdio>
#include "xvec.h"


// #define DEBUG  
#define MP

using namespace std;
// using namespace misc;
namespace PARAM{
    vector<string> pdblist;
    bool bypdb = false;

}


void print_help(int argc, char const *argv[]){
    std::cout << "       Get Protein or RNA Sequences from PDB, in Fasta Format" << endl;
    std::cout << "Usage: " << argv[0] <<" [options] pdb1[ pdb2 pdb3 ...]" << endl;
    std::cout << "                                                      " << endl;
    std::cout << "Options:" << endl;    
    std::cout << "             -bypdb     One sequence by pdb" << endl;
}

void read_params(int argc, char const *argv[]){
    if (argc < 2){
        print_help(argc, argv);
        exit(0);
    }

    PARAM::pdblist.clear();
    // for (int i0 = 1; i0 < argc ; i0++){
    //     PARAM::pdblist.push_back(argv[i0]);
    // }


    for(int i=1; i<argc; i++){
        string opt = argv[i];
        if (opt == "-bypdb"){
            PARAM::bypdb = true;
            continue;
        }
        else if (argv[i][0] == '-'){
            cerr << "unknown option: " << argv[i] << endl;
        }
        else {
            if (misc::file_exists(argv[i])){
                PARAM::pdblist.push_back(argv[i]);
                // cerr << "pdb added: " << argv[i]<< endl;
            }
            else {
                cerr << "File not exist: " << argv[i] << endl;
            }
            
        }
    }
}

double get_rmsd(rna p1, rna r2){


    return 0;

}

double get_rmsd(const structure & p1, const structure & p2){
    return 0;

}


double get_rmsd(std::vector<Xvec> v1, std::vector<Xvec> v2){

    if (v1.size() != v2.size()){
        cerr << "length not match: " << v1.size() << " " << v2.size() << endl;
        exit(1);
    }





    return 0;
}


int main(int argc, char const *argv[])
{

}
