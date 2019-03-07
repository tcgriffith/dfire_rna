// #include "PDB.h"

#include "dfire_PDB.h"

#include "dfire_dipolar.h"

#include "dfire_calculator.h"

// #include "dfire_rna.h"

// ========================
// main test
//
// ========================



// void test_smcraparam(){

//     smcra_param * param_p = new smcra_param();

//     cout << param_p->atomtype2[0] << endl;

//     string a = "A_C1'";
//     cout << param_p->atomtype[a] << endl;

// }

// void test_dihedral(){
//     string dfire_folder = "/home/tc/GIT/dfire_rna2/data/energyfiles";
//     dfire_dihedral ddih(dfire_folder);
//     string apdb = "/home/tc/GIT/dfire_rna2/test/1a9nR.pdb";
//     structure astr(apdb);
//     astr.rna_init_dihedrals();
//     double energy_dih = ddih.get_dfire(astr);
//     cout << energy_dih << endl;
// }

// void test_dih_init(){
//     string dfire_folder = "/home/tc/GIT/dfire_rna2/data/energyfiles";
//     dfire_dihedral ddih(dfire_folder);
// }

// void test_dipo_init(){
//     string dfire_folder = "/home/tc/GIT/dfire_rna2/data/energyfiles";
//     dfire_dipolar ddipo(dfire_folder);

//     ddipo.print_edfire_dipo(stdout, stdout);
// }

// void test_dipo(){
//     string dfire_folder = "/home/tc/GIT/dfire_rna2/data/energyfiles";
//     dfire_dipolar dc(dfire_folder);
//     string apdb = "/home/tc/GIT/dfire_rna2/test/1a9nR.pdb";
//     structure astr(apdb);
//     astr.rna_init_dihedrals();
//     astr.rna_init_polar();
//     std::pair<double, double> energy_dipo = dc.get_dfire(astr);

//     printf("%8.6f %8.6f\n", energy_dipo.first, energy_dipo.second);

// }


// void test_pvec(){
//     string apdb = "/home/tc/GIT/dfire_rna2/test/1a9nR.pdb";
//     structure astr(apdb);

//     astr.rna_init_polar();

//     for (chain & achain : astr.chains){
//         for (residue & ares : achain.residues){
//             for (atom & aatom : ares.atoms){
//                 if (aatom.is_polar()){
//                     printf("%-6s %-6s %8s\n",ares.get_resname().c_str(), aatom.get_name().c_str(),aatom.get_pvec().c_str());
//                 }
//             }
//         }
//     }



// }

// void test_idsd(){
//     string apdb = "/home/tc/GIT/dfire_rna2/test/1a9nR.pdb";
//     // string apdb = "/home/tc/GIT/spalign/SPalignNS/spalign_data/RMalign/datasets/PDB_3m/pdb/3g8t.pdb";
//     structure astr(apdb);

//     for(auto & achain :astr.chains){
//         // cout << achain.get_chaintype() <<endl;
//         // cout << achain.get_seq() << endl;
//         for (auto & ares: achain.residues){
//             cout << ares.get_resname() << endl;
//         }
//     }

// }

void test_pdb() {
    string apdb = "/home/tc/GIT/dfire_rna2/test/1a9nR.pdb";

    structure astr(apdb);

    for (auto &achain : astr.chains) {
        // cout << achain.get_chaintype() <<endl;
        // cout << achain.get_seq() << endl;
        for (auto &ares : achain.residues) {
            cout << ares.get_resname() << endl;
        }
    }


    // cout << "Hello world" << endl;
}

void test_dfirepdb() {

    string apdb = "/home/tc/GIT/dfire_rna2/test/1a9nR.pdb";

    // structure astr(apdb);
    rna astr(apdb);

    for (auto &achain : astr.chains) {
        // cout << achain.get_chaintype() <<endl;
        cout << achain.get_seq() << endl;
        // for (auto & ares: achain.residues){
        //     cout << ares.get_resname() << endl;
        // }
    }

}

void test_dih() {
    string apdb = "/home/tc/GIT/dfire_rna2/test/1a9nR.pdb";

    // structure astr(apdb);
    rna astr(apdb);

    astr.rna_init_dihedrals();

    for (auto &achain : astr.chains) {
        // cout << achain.get_chaintype() <<endl;
        // cout << achain.get_seq() << endl;
        for (auto &ares : achain.residues) {
            // cout << ares.get_resname() << endl;
            cout << ares.dihedrals_to_string() << endl;
        }
    }
}

void test_dipo() {
    string apdb = "example/1a9nR.pdb";
    rna astr(apdb);

    astr.rna_init_polar();

    for (auto &achain : astr.chains) {
        for (auto &ares : achain.residues) {
            if (ares.lack_atoms.size() > 0) {
                cout << ares.lack_atoms[0] << endl;
            }
            for (auto &aatom : ares.atoms) {
                // cout << aatom.get_name() << " " << aatom.is_polar() << endl;
                cout << aatom.get_name() << " " << aatom.get_pvec()  << endl;


            }
        }
    }

}

void test_polartypes() {
    for (auto &term : pdb_utils::atomtype) {
        if (pdb_utils::is_rna_polar_atom(term.first)) {
            cout << term.first << " " << term.second << endl;
        }
    }
}

void test_dfire_dipo() {
    string dfire_folder = "/home/tc/GIT/dfire_rna2/data/energyfiles";

    dfire_dipolar ddipo(dfire_folder);

    //

    cout << ddipo.dfire_dipo_pn.shape()[0] << endl;
    cout << ddipo.dfire_dipo_pn.shape()[1] << endl;
    cout << ddipo.dfire_dipo_pn.shape()[2] << endl;
    cout << ddipo.dfire_dipo_pn.shape()[3] << endl;

    string etable_pn = "/home/tc/GIT/dfire_rna2/test/trainxp/energy.dipo.pn2.dat";
    string etable_pp = "/home/tc/GIT/dfire_rna2/test/trainxp/energy.dipo.pp2.dat";



    // cout << ddipo.dfire_dipo_pn.shape()[0] << endl;

    ddipo.print_edfire_dipo(etable_pn, etable_pp);





}

void test_score_dipo() {
    string dfire_folder = "/home/tc/GIT/dfire_rna2/data/energyfiles";
    dfire_dipolar ddipo(dfire_folder);

    string apdb = "example/1a9nR.pdb";
    rna astr(apdb);

    astr.rna_init_polar();

    std::pair<double, double> score = ddipo.calc_dfire_dipo(astr);

    cout << score.first << " " << score.second << endl;


}

void train_dis() {
    // string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.l";
    // string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.tail5";
    string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.tail100";
    FILE *count_out = fopen("/home/tc/GIT/dfire_rna2/test/train_params/count.dat", "w");
    FILE *etable_out = fopen("/home/tc/GIT/dfire_rna2/test/train_params/energy.dat", "w");

    // dfire_dihedral * dc = new dfire_dihedral();
    dfire_calculator *dc = new dfire_calculator();


    double alpha = 1.61;
    double bin_size = 0.3;
    double max_dist = 15.0;

    char sout[200];

    sprintf(sout, "%.2f_%.2f_%.2f", alpha, bin_size, max_dist);

    // string myparam = "" + to_string(bin_size) +"_"+ to_string(max_dist);

    cout << sout << endl;



    dc->train_dfire(pdbl, alpha, bin_size, max_dist);

    dc->print_count(count_out);
    dc->print_edfire(etable_out);

    fclose(count_out);
    fclose(etable_out);

}

void test_score_dist() {

    string dfire_folder = "/home/tc/GIT/dfire_rna2/dataset/train_params/";
    // dfire_dipolar ddipo(dfire_folder);


    // always use pointer for dc

    dfire_calculator *dc = new dfire_calculator(dfire_folder);
    // dfire_calculator dc;
    // dc.init_dfire_etable(dfire_folder);
    // dfire_calculator * dc = new dfire_calculator();

    // dc->init_dfire_etable(dfire_folder);

    // dc->print_edfire(stdout);

    string apdb = "example/1a9nR.pdb";
    rna astr(apdb);

    cout << astr.get_natom() << endl;

    double score = dc->get_dfire(astr);

    cout << score << endl;

    // printf("%f\n", score);

    // // astr.rna_init_polar();

    // // std::pair<double, double> score = ddipo.calc_dfire_dipo(astr);

    // // cout << score.first << " " << score.second << endl;

}

void test_3b58abc() {
    string apdb = "example/3b58ABC.pdb";
    structure astr(apdb);

    // cout << astr.get_natom() << endl;

    // for (auto &achain : astr.chains) {

    //     // boost::filesystem::path p(apdb);

    //     string name = misc::basename(apdb);

    //     string fasta_info = "> " + name + " " + achain.get_chainid() + " " + to_string(achain.get_residue_num());

    //     cout << fasta_info << endl;
    //     cout << achain.get_seq() << endl;
    // }

}

void test_guess_ctype(){

    std::vector<std::string> tests = {
        "A",
        "VAL",
        "T",
        "5MT",
    };

    for (auto & a_test: tests){
        printf("%s %d\n", a_test.c_str(), pdb_utils::guess_ctype(a_test) );
    }



}

void test_folder_exist(){
    if (misc::file_exists("src/tester.cc")){
        cout << "yes" << endl;
    }
}


int main(int argc, char const *argv[]) {
    test_folder_exist();

    // test_guess_ctype();
    // test_3b58abc();
    // test_score_dist();

    // train_dis();
    // test_score_dipo();
    // test_dfire_dipo();
    // test_polartypes();
    // test_dipo();
    // test_dih();
    // test_dfirepdb();
    // test_idsd();
    // test_dipo();
    // test_pvec();

    // read_param(argc, argv);

    // cout << PARAM::dfire_folder << endl;

    // for (string & s : PARAM::pdblist){
    //     cout << s << endl;
    // }

    // run_dfire(argc, argv);
    /* code */
    return 0;
}

