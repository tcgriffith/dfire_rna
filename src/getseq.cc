// #include "structure.h"
// #include <boost/filesystem.hpp>

#include "PDB.h"
#include <algorithm>
#include <unordered_set>
#include <cstdio>


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

void getseq_by_chain(){

    int ne = PARAM::pdblist.size();
    vector<string> sdim (ne, "");
    FILE *fp = stdout;

#if defined MP
    #pragma omp parallel for
#endif
    for (int i = 0; i < ne; i++){
        string apdb =PARAM::pdblist[i];
        structure * astr = new structure(apdb);
        char sout[8000] = "";
        for (auto & achain: astr->chains){

            // boost::filesystem::path p(apdb);

            string name = misc::basename(apdb);

            // char tmp[8000];

            sprintf(sout, "%s> %s %s %d\n%s\n", sout,name.c_str(), achain.get_chainid().c_str(), achain.get_residue_num(),achain.get_seq().c_str());

            // string fasta_info = "> " + name + " "+achain.get_chainid() + " " +to_string(achain.get_residue_num());

            // cout << fasta_info << endl;
            // cout << achain.get_seq() << endl;
        }
        #if defined MP
        sdim[i] = sout;

        #else
        fprintf(fp,"%s",sout);
        #endif

        delete astr;
    }

#if defined MP
    for (int i = 0 ; i < ne; i++){
        fprintf(fp, "%s", sdim[i].c_str());
    }
#endif    

}

void getseq_by_pdb(){

    

    for (string apdb : PARAM::pdblist){
        std::unordered_set<string> ctypes = {};
        structure * astr = new structure(apdb);

        string name = misc::basename(apdb);
        string fasta_info = "> " + name + " " +to_string(astr->get_nres());
        string seq = "";

        for (auto & achain: astr->chains){
            // boost::filesystem::path p(apdb);
            if (!ctypes.count(achain.get_chaintype())){
                ctypes.insert(achain.get_chaintype());
            }
            seq = seq + achain.get_seq();


        }

        cout << fasta_info << endl;
        cout << seq << endl;
        // cout << ctypes.size() << endl;

        if (ctypes.size() > 1){
            cerr << "# Warning: mixed chain types" << endl;
        }
        
        delete astr;
    }   



}

void run_getseq(int argc, char const *argv[]){

    read_params(argc, argv);

    // cout << PARAM::bypdb << endl;

    if (PARAM::bypdb) {
        getseq_by_pdb();
    }
    else {
        getseq_by_chain();
    }
}

int main(int argc, char const *argv[])
{
    run_getseq(argc,argv);
    /* code */
    return 0;
}
