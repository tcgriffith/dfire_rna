#include "dfire_dihedral.h"

void dfire_dihedral::print_edfire_dih(FILE * out) {
    // FILE * out;
    // out = fopen(out.c_str(),"w");

    for (int i=0;i<dtype;i++){
        fprintf(out,"type_%d",i);
        for(int j=0;j<dbin;j++){
                fprintf(out,"\t%8.6f",dfire_dih[i][j]);
        }
        fprintf(out,"\n");
    }
 
}


void dfire_dihedral::print_count_dih(FILE * out){
    for (int i=0;i<dtype;i++){
        fprintf(out,"type_%d",i);
        for(int j=0;j<dbin;j++){
                fprintf(out,"\t%d",count_dih[i][j]);
        }
        fprintf(out,"\n");
    }
}


void dfire_dihedral::train_dfire_dih(string pdbl_fpath, int length_limit){


    string train_list=pdbl_fpath;

    std::ifstream fs(train_list.c_str());
    std::string line;

    while(std::getline(fs, line)){

        cout << line <<endl;
        // structure *astr = new structure();
        rna *astr = new rna(); 
        astr->readpdb(line);
        cout << astr->get_nres() <<endl;

        // remove > 500 nt structures

        if (astr->get_nres() > length_limit){
            cout << "skipped" << endl;
            continue;
        }
        astr->rna_init_dihedrals(); 

        count_dihs(*astr);
    }

    count_to_dfiredih();

}

void dfire_dihedral::init_dfire_etable(string energy_dir){

    std::string dfire_file = energy_dir + "/energy.dih.dat";

    // if (! boost::filesystem::exists(dfire_file)){
    if (! misc::file_exists(dfire_file)){
        cerr << "# File not found:" << dfire_file.c_str() << endl;
        exit(1);
    }

    std::ifstream fs(dfire_file);
    std::string line;
    string ss[2];

    int row=0;

    while (std::getline(fs, line)){
        std::istringstream iss(line);
        iss >> ss[0];
        for (int m = 0; m < dbin; m++)
        {
            iss >> dfire_dih[row][m];
        }
        row++;
    }
};

double dfire_dihedral::calc_dfire_dih(rna &astr){ 

    double total_energy=0.0;

    for (auto & achain : astr.chains){
        if (achain.get_chaintype() != "NT") continue;
        for (auto & ares : achain.residues){ 
            if (ares.has_dih()){
                double *dihs = ares.get_dihedrals(); 
                for (int i = 0; i < 7; i++){
                    double tmpdih=dihs[i];
                    if (tmpdih < 0) tmpdih = 360 + tmpdih;
                    int dihbin = (int)tmpdih/this->dbin_width;
                    total_energy += dfire_dih[i][dihbin];
                }
            }
        }
    }

    return total_energy;
}

// double dfire_dihedral::dihedral(atom &aa, atom &ab, atom &ac, atom &ad){
//     Xvec x1(aa.get_x());
//     Xvec x2(ab.get_x());
//     Xvec x3(ac.get_x());
//     Xvec x4(ad.get_x());

//     double tor = torsion_xvec(x1,x2,x3,x4);

//     return tor;

// }

void dfire_dihedral::count_dihs(rna & astr){ 

    for (auto & achain : astr.chains){
        for (auto & ares : achain.residues){ 
            if (ares.has_dih()){
                double *dihs = ares.get_dihedrals();
                for (int i = 0; i < 7; i++){
                    double tmpdih=dihs[i];
                    if (tmpdih < 0) tmpdih = 360 + tmpdih;

                    int dihbin = (int)tmpdih/this->dbin_width;

                    count_dih[i][dihbin]++;
                }
            }
        }
    }
}



void dfire_dihedral::count_to_dfiredih(){

    double dih_colsum[dbin]= {};
    double dih_rowsum[dtype]={};
    double dih_total=0;

    for (int i = 0; i < dtype; i++){
        for (int j = 0; j < dbin; j++){
            dih_colsum[j] += count_dih[i][j];
            dih_rowsum[i] += count_dih[i][j];
            dih_total += count_dih[i][j];
        }
    }

    if (this->dih_ref_type == 0) {
        for (int i = 0; i < dtype; i++){
            for (int j = 0; j < dbin; j++){
                double quot, out;

                quot = 1.0 * count_dih[i][j] * dih_total / dih_rowsum[i] / dih_colsum[j];

                if (quot < 0.000001) out = 10.0;
                else out = - 0.001987 * 300 * log(quot);

                dfire_dih[i][j]=out;
            }
        }
    }
    else if (this->dih_ref_type == 1){
        for (int i = 0; i < dtype; i++){
            for (int j = 0; j < dbin; j++){
                double quot, out;

                quot = 1.0 * count_dih[i][j] / dih_rowsum[i] * 1.0 / dbin;

                if (quot < 0.000001) out = 10.0;
                else out = - 0.001987 * 300 * log(quot);

                dfire_dih[i][j]=out;
            }
        }
    }

}

