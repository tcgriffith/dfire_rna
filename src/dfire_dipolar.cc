#include "dfire_dipolar.h"


// todo : get_dfire()

// #define DEBUG



std::pair<double, double> dfire_dipolar::calc_dfire_dipo(rna &astr){

    double etotal_pn = 0.;

    double etotal_pp = 0.;

    for (auto & achain : astr.chains){
        if (achain.get_chaintype() != "NT") continue;
        for (auto & bchain : astr.chains){
            if (bchain.get_chaintype() != "NT") continue;
            if (achain.get_chainid() > bchain.get_chainid()) continue;
            for (auto & ares : achain.residues){
                if (ares.lack_atoms.size() > 0) continue;
                for (auto & bres: bchain.residues){    
                    if (bres.lack_atoms.size() > 0) continue;
                    if (achain.get_chainid() == bchain.get_chainid()){
                        if (ares.get_residsd() >= bres.get_residsd()) continue;
                    }
                    for(auto & aatom : ares.atoms){
                        // atom_num++;
                        for (auto & batom: bres.atoms){

                            if ((!aatom.is_polar()) && (!batom.is_polar())) continue;

                            double d = aatom.distance(batom);

                            int bin = int(d * 2);
                            if (bin >= 30) continue;

                            int atype_ia=-1, atype_ib=-1;

                            if (aatom.is_polar()){

                                atype_ia = pdb_utils::get_polartype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_atomtype_int(batom.get_type());
                                if (atype_ia == -1 || atype_ib == -1) continue;

                                double cos_a =aatom.get_cos_angle_p(batom);

                                int abin = cosang2bin(cos_a);
                                if (abs(cos_a) > 1)  {
                                    cerr << ares.get_residsd() << " " << aatom.get_name()<< bres.get_residsd() << " " << batom.get_name() << endl;
                                    exit(1);
                                }

                                // count_dipo_pn[atype_ia][atype_ib][abin][bin]++;


                                #if defined DEBUG

                                printf("pn %d %s %d %s %8.4f\n", ares.get_residsd(), aatom.get_name().c_str(), bres.get_residsd(), batom.get_name().c_str(), dfire_dipo_pn[atype_ia][atype_ib][abin][bin]);

                                #endif

                                etotal_pn = etotal_pn + dfire_dipo_pn[atype_ia][atype_ib][abin][bin];

                            }

                            if (batom.is_polar()){

                                atype_ia = pdb_utils::get_atomtype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_polartype_int(batom.get_type());
                                if (atype_ia == -1 || atype_ib == -1) continue;

                                double cos_b =batom.get_cos_angle_p(aatom);

                                if (abs(cos_b) > 1)  {
                                    cerr << ares.get_residsd() << " " << aatom.get_name()<< bres.get_residsd() << " " << batom.get_name() << endl;
                                    exit(1);
                                }

                                int abin = cosang2bin(cos_b);
                                #if defined DEBUG

                                printf("np %d %s %d %s %8.4f\n", ares.get_residsd(), aatom.get_name().c_str(), bres.get_residsd(), batom.get_name().c_str(), dfire_dipo_pn[atype_ib][atype_ia][abin][bin]);

                                #endif

                                // count_dipo_pn[atype_ib][atype_ia][abin][bin]++;
                                etotal_pn = etotal_pn + dfire_dipo_pn[atype_ib][atype_ia][abin][bin];
                            }

                            if (aatom.is_polar() && batom.is_polar()){

                                atype_ia = pdb_utils::get_polartype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_polartype_int(batom.get_type());
                                if (atype_ia == -1 || atype_ib == -1) continue;


                                double cos_ab = aatom.get_cos_angle_pq(batom);

                                int abbin =cosang2bin(cos_ab);
                                if (abs(cos_ab) > 1)  {
                                    cerr << ares.get_residsd() << " " << aatom.get_name()<< bres.get_residsd() << " " << batom.get_name() << endl;
                                    exit(1);
                                }

                                #if defined DEBUG

                                printf("pp %d %s %d %s %8.4f\n", ares.get_residsd(), aatom.get_name().c_str(), bres.get_residsd(), batom.get_name().c_str(), dfire_dipo_pp[atype_ia][atype_ib][abbin][bin]);

                                #endif

                                etotal_pp = etotal_pp + dfire_dipo_pp[atype_ia][atype_ib][abbin][bin];



                                // count_dipo_pp[atype_ia][atype_ib][abbin][bin]++;
                                // count_dipo_pp[atype_ib][atype_ia][abbin][bin]++;

                            }

                        }
                    }

                }
            }
        }
    }

    std::pair<double, double> etotal(etotal_pn,etotal_pp);

    return etotal;

}

void dfire_dipolar::init_dfire_etable_dipo(string energy_dir){
    std::string dfire_file1 = energy_dir + "/energy.dipo.pn.dat";
    std::string dfire_file2 = energy_dir + "/energy.dipo.pp.dat";

    if (! misc::file_exists(dfire_file1)){
        cerr << "# File not found: " << dfire_file1.c_str() << endl;
        exit(1);
    }
    if (! misc::file_exists(dfire_file2)){
        cerr << "# File not found: " << dfire_file2.c_str() << endl;
        exit(1);
    }

    std::ifstream fs(dfire_file1);
    std::string line;
    string ss[3];

    int row=0;
// pn

    while (std::getline(fs, line)){
        std::istringstream iss(line);
        iss >> ss[0] >> ss[1] >> ss[2];

        int id1 = pdb_utils::get_polartype_int(ss[0]);
        int id2 = pdb_utils::get_atomtype_int(ss[1]);
        if (id1 < 0) continue; // atom1 is not polar atom
        int id_dipo = atof(ss[2].c_str());

        for (int id_mbin = 0; id_mbin < mbin; id_mbin ++){
            iss >> dfire_dipo_pn[id1][id2][id_dipo][id_mbin];
        }
    }

// pp
    fs.close();

    std::ifstream fs2(dfire_file2);
    while (std::getline(fs2, line)){
        std::istringstream iss(line);
        iss >> ss[0] >> ss[1] >> ss[2];

        int id1 = pdb_utils::get_polartype_int(ss[0]);
        if (id1 < 0) continue; // atom1 is not polar atom
        int id2 = pdb_utils::get_polartype_int(ss[1]);
        if (id2 < 0) continue; // atom2 is not polar atom

        int id_dipo = atof(ss[2].c_str());

        for (int id_mbin = 0; id_mbin < mbin; id_mbin ++){
            iss >> dfire_dipo_pp[id1][id2][id_dipo][id_mbin];
        }
    }

    fs2.close();

}





void dfire_dipolar::init_dipo_arrays(){
    for(index_i i = 0; i < polartypen; i ++)
        for (index_i j = 0; j < matype; j ++)
            for (index_i k = 0; k < dipo_bin; k ++)
                for (index_i l = 0; l < mbin;l ++){
                    count_dipo_pn[i][j][k][l] = 0;
                    dfire_dipo_pn[i][j][k][l] = 0;
                }

    for(index_d i = 0; i < polartypen; i ++)
        for (index_d j = 0; j < polartypen; j ++)
            for (index_d k = 0; k < dipo_bin; k ++)
                for (index_d l = 0; l < mbin;l ++){
                    count_dipo_pp[i][j][k][l] = 0;
                    dfire_dipo_pp[i][j][k][l] = 0;
                }
}

void dfire_dipolar::train_dfire_dipo(std::string pdbl_fpath, int length_limit){
    this->init_dipo_arrays();
    std::string train_list=pdbl_fpath;

    std::ifstream fs(train_list.c_str());
    std::string line;

    while(std::getline(fs, line)){

        cout << line <<endl;
        rna *astr = new rna();

        astr->readpdb(line);
        cout << astr->get_nres() <<endl;
        // astr->chains[0].rna_init_polar();
        if (astr->get_nres() > length_limit){
            cout << "skipped" << endl;
            continue;

        }
        if (astr->rna_init_polar()){
            cout << "init polar failed, skipped" << endl;
            continue;

        }
        count_dipo(*astr);
        // count_dihs(astr->chains[0]);
    }

    count_to_dfiredipo();
}




void dfire_dipolar::count_dipo(rna &astr){

    for (auto & achain : astr.chains){
        for (auto & bchain: astr.chains){
            for (auto & ares : achain.residues){
                if (ares.lack_atoms.size() > 0) continue;
                for (auto & bres: bchain.residues){    
                    
                    if (bres.lack_atoms.size() > 0) continue;
                    if (achain.get_chainid() == bchain.get_chainid()){
                        if (ares.get_residsd() >= bres.get_residsd()) continue;
                    }

                    if (achain.get_chainid() > bchain.get_chainid()) continue;

                    for(auto & aatom : ares.atoms){
                        // atom_num++;
                        for (auto & batom: bres.atoms){

                            if (!aatom.is_polar() && !batom.is_polar()) continue;

                            double d2 = aatom.distance2(batom);
                            if (d2 > 225) continue;

                            double d = aatom.distance(batom);

                            int bin = int(d * 2);
                            if (bin >= 30) continue;

                            int atype_ia=-1, atype_ib=-1;

                            // atype_ia = pdb_utils::get_atomtype_int(aatom.get_type());
                            // atype_ib = pdb_utils::get_atomtype_int(batom.get_type());

                            // if (atype_ia == -1 || atype_ib == -1) continue;


                            if (aatom.is_polar()){
                                atype_ia = pdb_utils::get_polartype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_atomtype_int(batom.get_type());
                                if (atype_ia == -1 || atype_ib == -1) continue;


                                double cos_a =aatom.get_cos_angle_p(batom);


                                int abin = cosang2bin(cos_a);

                                count_dipo_pn[atype_ia][atype_ib][abin][bin]++;

                            }

                            if (batom.is_polar()){

                                atype_ia = pdb_utils::get_atomtype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_polartype_int(batom.get_type());
                                if (atype_ia == -1 || atype_ib == -1) continue;

                                double cos_b =batom.get_cos_angle_p(aatom);

                                int abin = cosang2bin(cos_b);

                                count_dipo_pn[atype_ib][atype_ia][abin][bin]++;
                            }

                            if (aatom.is_polar() && batom.is_polar()){

                                atype_ia = pdb_utils::get_polartype_int(aatom.get_type());
                                atype_ib = pdb_utils::get_polartype_int(batom.get_type());
                                if (atype_ia == -1 || atype_ib == -1) continue;

                                double cos_ab = aatom.get_cos_angle_pq(batom);

                                int abbin =cosang2bin(cos_ab);

                                count_dipo_pp[atype_ia][atype_ib][abbin][bin]++;
                                count_dipo_pp[atype_ib][atype_ia][abbin][bin]++;

                            }

                            // std::cout << aatom.get_type() << " " <<aatom.get_type() <<endl;


                        }
                    }


                }
            }
        }
    }


}

// angle in radian
//




void dfire_dipolar::count_to_dfiredipo(){

    for (index_d i = 0; i < count_dipo_pn.shape()[0]; i++){
        for(index_d j = 0; j < count_dipo_pn.shape()[1]; j++){
            for (index_d k = 0; k < count_dipo_pn.shape()[2]; k++){
                for (index_i l = 0; l < count_dipo_pn.shape()[3]; l++){
                    double quot2, out2, r_rat;
                    // double quot, , out;
                    r_rat = pow((l/2+0.25)/14.75, ALPHA);

                    if (r_rat < 1e-6) quot2 = 0.0;
                    else quot2 = count_dipo_pn[i][j][k][l]/(r_rat *count_dipo_pn[i][j][k][29]);

                    if (quot2 < 1e-6) out2 = 10.0;
                    else out2 =  -.001987 *300 * log(quot2);
                    dfire_dipo_pn[i][j][k][l] = out2;
                }
            }
        }
    }



    for (index_i i = 0; i < count_dipo_pp.shape()[0]; i++){
        for(index_i j = 0; j < count_dipo_pp.shape()[1]; j++){
            for (index_i k = 0; k < count_dipo_pp.shape()[2]; k++){
                for (index_i l = 0; l < count_dipo_pp.shape()[3]; l++){

                    //  dipolar - dipolar
                    double quot, r_rat, out;
                    r_rat = pow((l/2+0.25)/14.75, ALPHA);

                    if (r_rat < 1e-6) quot = 0.0;
                    else quot = count_dipo_pp[i][j][k][l]/(r_rat *count_dipo_pp[i][j][k][29]);

                    if (quot < 1e-6) out = 10.0;
                    else out =  -.001987 *300 * log(quot);
                    dfire_dipo_pp[i][j][k][l] = out;
                }
            }

        }
    }
}

void dfire_dipolar::print_count_dipo(FILE * out_pn, FILE * out_pp){

// pn
    for (index_d i = 0; i < count_dipo_pn.shape()[0]; i++){
        for(index_d j = 0; j < count_dipo_pn.shape()[1]; j++){
            for (index_d k = 0; k < count_dipo_pn.shape()[2]; k++){
                fprintf(out_pn,"%s\t%s\t%d",pdb_utils::typepolar[i].c_str(),pdb_utils::typeatom[j].c_str(),int(k));
                for (index_i l = 0; l < count_dipo_pn.shape()[3]; l++){
                    fprintf(out_pn,"\t%d",count_dipo_pn[i][j][k][l]);
                }
                fprintf(out_pn,"\n");
            }

        }
    }
// pp
    for (index_d i = 0; i < count_dipo_pp.shape()[0]; i++){
        for(index_d j = 0; j < count_dipo_pp.shape()[1]; j++){
            for (index_d k = 0; k < count_dipo_pp.shape()[2]; k++){
                fprintf(out_pp,"%s\t%s\t%d",pdb_utils::typepolar[i].c_str(),pdb_utils::typepolar[j].c_str(),int(k));
                for (index_i l = 0; l < count_dipo_pp.shape()[3]; l++){
                    fprintf(out_pp,"\t%d",count_dipo_pp[i][j][k][l]);
                }
                fprintf(out_pp,"\n");
            }

        }
    }

}

void dfire_dipolar::print_count_dipo(std::string & sout_pn, std::string & sout_pp){
    FILE * out_pn;
    out_pn = fopen(sout_pn.c_str(),"w");

    FILE * out_pp;
    out_pp = fopen(sout_pp.c_str(),"w");

    print_count_dipo(out_pn, out_pp);

    fclose(out_pn);
    fclose(out_pp);
}

void dfire_dipolar::print_edfire_dipo(std::string & sout_pn, std::string & sout_pp){
    FILE * out_pn;
    out_pn = fopen(sout_pn.c_str(),"w");

    FILE * out_pp;
    out_pp = fopen(sout_pp.c_str(),"w");

    print_edfire_dipo(out_pn, out_pp);

    fclose(out_pn);
    fclose(out_pp);
}

void dfire_dipolar::print_edfire_dipo(FILE * out_pn, FILE * out_pp){

    for (index_d i = 0; i < dfire_dipo_pn.shape()[0]; i++){
        for(index_d j = 0; j < dfire_dipo_pn.shape()[1]; j++){
            for (index_d k = 0; k < dfire_dipo_pn.shape()[2]; k++){
                fprintf(out_pn,"%s\t%s\t%d",pdb_utils::typepolar[i].c_str(),pdb_utils::typeatom[j].c_str(),int(k));
                for (index_d l = 0; l < dfire_dipo_pn.shape()[3]; l++){
                    fprintf(out_pn,"\t%8.6f",dfire_dipo_pn[i][j][k][l]);
                }
                fprintf(out_pn,"\n");
            }

        }
    }


    for (index_d i = 0; i < dfire_dipo_pp.shape()[0]; i++){
        for(index_d j = 0; j < dfire_dipo_pp.shape()[1]; j++){
            for (index_d k = 0; k < dfire_dipo_pp.shape()[2]; k++){
                fprintf(out_pp,"%s\t%s\t%d",pdb_utils::typepolar[i].c_str(),pdb_utils::typepolar[j].c_str(),int(k));
                for (index_d l = 0; l < dfire_dipo_pp.shape()[3]; l++){
                    fprintf(out_pp,"\t%8.6f",dfire_dipo_pp[i][j][k][l]);
                }
                fprintf(out_pp,"\n");
            }

        }
    }

}



int dfire_dipolar::ang2bin(double angle){
    int angle_bin=-1;
    double cos_angle = cos(angle);
    angle_bin = (int) (cos_angle + 1.0)/ (1.0/3);
    return angle_bin;
}


// angle bins, cos(angle)

// -1   ~ -2/3 -> 0
// -2/3 ~ -1/3 -> 1
// -1/3 ~  0   -> 2
// 0    ~  1/3 -> 3
// 1/3  ~  2/3 -> 4
// 2/3  ~  1   -> 5



int dfire_dipolar::cosang2bin(double cos_angle){
    int angle_bin=-1;
    angle_bin = (int) ((cos_angle + 1.0)/ (1.0/3));
    if (angle_bin == 6){
        angle_bin = 5;
    }
    return angle_bin;
}

