// #include "structure.h"
// #include <boost/filesystem.hpp>

#include "PDB.h"

using namespace std;
// using namespace misc;

// add -merge option, output a unique sequence for a pdb

namespace params{
    vector<string> pdblist;
    bool merge = false;
}


void read_params(int argc, char const *argv[]){
    if (argc < 2){
        std::cout << "       Get Protein or RNA Sequences from PDB, in Fasta Format" << endl;
        std::cout << "Usage: " << argv[0] <<" pdb1[ pdb2 pdb3 ...]" << endl;
        exit(0);
    }



    params::pdblist.clear();
    for (int i0 = 1; i0 < argc ; i0++){
        string opt = argv[i0];
        if (opt == "-merge") {
            params::merge = true;
            // i0 ++;
            continue;
        }
        params::pdblist.push_back(argv[i0]);
    }
}

void run_getseq(int argc, char const *argv[]){

    read_params(argc, argv);

    if (params::merge){
        for (string apdb : params::pdblist){
            structure * astr = new structure(apdb);
            string seq="";

            string name = misc::basename(apdb);

            string fasta_info= "> " + name + " "+ to_string(astr->get_nres()) + " ";

            for (auto & achain: astr->chains){
                seq = seq + achain.get_seq();
                fasta_info = fasta_info + achain.get_chainid() + " ";
            }

            cout << fasta_info << endl;
            cout << seq << endl;
        }
    }
    else {
        for (string apdb : params::pdblist){
            structure * astr = new structure(apdb);
            for (auto & achain: astr->chains){

                // boost::filesystem::path p(apdb);

                string name = misc::basename(apdb);

                string fasta_info = "> " + name + " "+achain.get_chainid() + " " +to_string(achain.get_residue_num());

                cout << fasta_info << endl;
                cout << achain.get_seq() << endl;
            }
        }

    }




}

int main(int argc, char const *argv[])
{
    run_getseq(argc,argv);
    /* code */
    return 0;
}
