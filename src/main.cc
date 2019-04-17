#include "main.h"

// ========================
// main test
//
// ========================

#define MP

namespace PARAM{
    vector<string> pdblist;
    string dfire_folder;
    bool norm = false;
}
 

bool calc_dfire_list(vector<string> pdblist, string dfire_folder){

    // dfire_calculator dc(dfire_folder);
    dfire_calculator * dc = new dfire_calculator(dfire_folder);
    // dfire_dihedral ddih(dfire_folder); 
    // dfire_dipolar ddipo(dfire_folder);
    // string dfire_dih_folder = "/home/tc/GIT/dfire_rna2/tmpdata";
    // dc.init_dfire_dih(dfire_dih_folder);

    int ne = pdblist.size();
    vector<string> sdim (ne, "");
    FILE *fp = stdout;

#if defined MP
    #pragma omp parallel for
#endif
    for (int i = 0; i < ne; i++){
        string line = pdblist[i];
        // boost::filesystem::path p(line);



        // if (!boost::filesystem::exists(p)){
        if (!misc::file_exists(line)){
            cerr << "pdb " << line << " not exist, skipped" << endl;
            continue;
        }

        // structure *astr = new structure();
        rna *astr = new rna(line);
        // astr->readpdb(line);
        // astr->chains[0].init_dihedrals();
        // astr->rna_init_dihedrals();
        // astr->rna_init_polar();

        // double energy_td = dc.calc_3drna(*astr);
       

        double energy = dc->get_dfire(*astr);

        double energy_dipo_pn = 0.;
        double energy_dipo_pp = 0.;
        double energy_dih = 0.;

        {
        // double energy_dih = ddih.get_dfire(*astr); 
        // std::pair<double, double> energy_dipo = ddipo.get_dfire(*astr);
        // double energy_dipo_pn = 0.energy_dipo.first;
        // double energy_dipo_pp = 0.energy_dipo.second;
        }






        char sout[200];
        // sprintf(char *__restrict __s, const char *__restrict __format, ...)
        // sprintf(sout, "%s %8.6f %8.6f %8.6f %8.6f %8.6f \n", line.c_str(), energy, energy_dih, energy_dipo_pn, energy_dipo_pp, energy_td);
            // sprintf(sout, "%s %8.6f %8.6f %8.6f %8.6f \n", line.c_str(), energy, energy_dih, 0., 0.);

        sprintf (sout, "%s %8.6f ", line.c_str(), energy);
        // sprintf(sout, "%s %8.6f %8.6f %8.6f %8.6f ", line.c_str(), energy, energy_dih, energy_dipo_pn, energy_dipo_pp);
        // sprintf(sout, "%s %8.6f %8.6f %8.6f %8.6f \n", line.c_str(), energy, 0., 0., 0.);

        if (PARAM::norm){ // energy per residue

            double energy_norm_by_res = 1.0 * energy / astr->get_nres();
            double energy_norm_by_pair = 1.0 * energy / astr->n_atompairs;

            sprintf(sout, "%s %8.6f %8.6f \n", sout, energy_norm_by_res, energy_norm_by_pair);

            // double energy_dih_n = 1.0 * energy_dih / astr->get_nres();
            // double energy_dipo_pn_n = 1.0 * energy_dipo_pn / astr->get_nres();
            // double energy_dipo_pp_n = 1.0 * energy_dipo_pp / astr->get_nres();

            // sprintf(sout, "%s %8.6f %8.6f %8.6f %8.6f \n",sout, energy_norm, energy_dih_n, energy_dipo_pn_n, energy_dipo_pp_n);

        }
        else {
            sprintf(sout,"%s \n", sout);
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

    delete dc;

    return 0;
}

void print_help(char const *argv[]){
        std::cout << "#######################################################" << endl;
        std::cout << "# Calculate dfire_rna score for a pdb or a list of pdbs" << endl;
        std::cout << "#######################################################" << endl;
        std::cout << "Usage: " << argv[0] <<" pdb " << endl;
        std::cout << "   or: " << argv[0] <<" [ options ] " << endl;


        std::cout << "Options:" << endl;
        std::cout << "   pdb [ pdb2 pdb3 ...], input RNA structures in pdb format" << endl;
        std::cout << "   -d directory,         OPTIONAL, override default directory of energyfiles" << endl;
        std::cout << "                         default:" << PARAM::dfire_folder  << endl;
        std::cout << "   -norm                 normalize DFIRE score by RNA length (experimental)" << endl;
        std::cout << "   -l pdblist,           A list of absolute paths to pdb files (plain text) UTF-8 encoding" << endl;

        exit(1);
}

void read_param(int argc, char const *argv[]){

    string projfolder = getenv("DFIRE_RNA_HOME");

    PARAM::dfire_folder = projfolder+ "/data/energyfiles";

    // vector<string> pdblist;

    for(int i=1; i<argc; i++){
        string opt = argv[i];
        if (opt == "-d") {
            i++;
            PARAM::dfire_folder = argv[i];
        }
        else if (opt == "-l") {
            i++;

            std::ifstream fs(argv[i]);
            std::string line;

            while(std::getline(fs, line)){
                PARAM::pdblist.push_back(line);
            }
        }
        else if (opt == "-s"){
            continue;
        }
        else if (opt == "-norm"){
            PARAM::norm = true;
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

    if (argc < 2){
        print_help(argv);

    }
    if (PARAM::dfire_folder.length() < 3){
        cerr << "!! dfire_rna root not found: " << PARAM::dfire_folder  << endl;
    }
}



int run_dfire(int argc, char const *argv[])
{
    //#######################
    // find executable
    // https://stackoverflow.com/a/49227138/7503413
    //#######################

    read_param(argc, argv);

    return calc_dfire_list(PARAM::pdblist, PARAM::dfire_folder);
}


int main(int argc, char const *argv[])
{
    run_dfire(argc, argv);
    return 0;
}
