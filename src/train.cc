// #include "dfire_rna.h"

// #include "dfire_dihedral.h"

#include "dfire_calculator.h"

// #include "dfire_dipolar.h"

using namespace std;



void train_3drna(){

    // string pdbl = "/home/tc/GIT/dfire_rna2/test/pdb.l";
    string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.l";
    // string etable_out = "/home/tc/GIT/dfire_rna2/test/trainxp/energy.dih.dat";
    // string count_out = "/home/tc/GIT/dfire_rna2/test/trainxp/count.dih.dat";

    FILE * count_out =fopen("/home/tc/GIT/dfire_rna2/test/train_3drna/count.dat","w");
    FILE * etable_out = fopen("/home/tc/GIT/dfire_rna2/test/train_3drna/energy.dat","w");

    dfire_calculator * dc = new dfire_calculator();

    dc->train_3drna(pdbl);

    dc->print_count(count_out);
    dc->print_edfire(etable_out);

    fclose(count_out);
    fclose(etable_out);

}


void train_dis_181217(){
    string pdbl = "dataset/PDB/RNA_hires_xp/RNA_hires_xp_pdbs.trainpdbs";

    vector<double> valpha = {1.30, 1.50, 1.61, 1.70, 1.80, 2.0};
    // vector<double> valpha = {1.61};
    // vector<double> vbin_size = {0.5};
    vector<double> vbin_size = {0.3, 0.5, 0.7, 1.0, 1.5, 2.0 };
    vector<double> vmax_dist ={14.0, 15.0,16.0,17.0,18.0,19.0,20.0};
    // vector<double> vmax_dist ={18.0};



    for (auto& alpha: valpha){
        for (auto& bin_size: vbin_size){
            for (auto& max_dist: vmax_dist){

                char path0[200];
                sprintf(path0,"test/train_181217/%.2f_%.2f_%.2f", alpha, bin_size,max_dist);



                string command ="mkdir -p "+ string(path0); 
                system(command.c_str());

                string path1 = string(path0) + "/count.dat";
                string path2 = string(path0) + "/energy.dat";

                if (misc::file_exists(path1)){
                    continue;
                }

                FILE * count_out =fopen(path1.c_str(),"w");
                FILE * etable_out = fopen(path2.c_str(),"w");

                dfire_calculator * dc = new dfire_calculator();

                dc->train_dfire(pdbl, alpha, bin_size, max_dist);

                dc->print_count(count_out);
                dc->print_edfire(etable_out);

                fclose(count_out);
                fclose(etable_out);

                delete dc;

                cout << path0 << " finished!" << endl;

            }
        }
    }

    // vector<double> vmax_dist ={17.0,18.0,19.0,20.0};    

}



// alpha 1.61, 
// bin 0.3 0.5, 1.0, 
// max dist 15, 18, 20

void train_dis_181219(){
    string pdbl = "dataset/PDB/RNA_hires_xp/RNA_hires_xp_pdbs.trainpdbs";

    // vector<double> valpha = {1.30, 1.50, 1.61, 1.70, 1.80, 2.0};
    vector<double> valpha = {1.61};
    vector<double> vbin_size = {0.3, 0.5, 1.0};
    // vector<double> vbin_size = {0.3, 0.5, 0.7, 1.0, 1.5, 2.0 };
    // vector<double> vmax_dist ={14.0, 15.0,16.0,17.0,18.0,19.0,20.0};
    vector<double> vmax_dist ={15.0, 18.0, 20};



    for (auto& alpha: valpha){
        for (auto& bin_size: vbin_size){
            for (auto& max_dist: vmax_dist){

                char path0[200];
                sprintf(path0,"test/train_181219/%.2f_%.2f_%.2f", alpha, bin_size,max_dist);



                string command ="mkdir -p "+ string(path0); 
                system(command.c_str());

                string path1 = string(path0) + "/count.dat";
                string path2 = string(path0) + "/energy.dat";

                if (misc::file_exists(path1)){
                    continue;
                }

                FILE * count_out =fopen(path1.c_str(),"w");
                FILE * etable_out = fopen(path2.c_str(),"w");

                dfire_calculator * dc = new dfire_calculator();

                dc->train_dfire(pdbl, alpha, bin_size, max_dist);

                dc->print_count(count_out);
                dc->print_edfire(etable_out);

                fclose(count_out);
                fclose(etable_out);

                delete dc;

                cout << path0 << " finished!" << endl;

            }
        }
    }

    // vector<double> vmax_dist ={17.0,18.0,19.0,20.0};    

}


void train_dis_190102(){
    string pdbl = "dataset/PDB/RNA_hires_xp/RNA_hires_xp_pdbs.trainpdbs";

    // vector<double> valpha = {1.30, 1.50, 1.61, 1.70, 1.80, 2.0};
    vector<double> valpha = {1.61};
    // vector<double> vbin_size = {0.3, 0.5, 1.0};
    // vector<double> vbin_size = {0.5};
    vector<double> vbin_size = {0.3, 0.5, 0.7, 0.9};
    vector<double> vmax_dist ={13.0, 15.0, 17.0, 18.0,  19.0, 21.0};
    // vector<double> vmax_dist ={15.0, 18.0, 20};
    // vector<double> vmax_dist ={15.0};


    for (auto& alpha: valpha){
        for (auto& bin_size: vbin_size){
            for (auto& max_dist: vmax_dist){

                char path0[200];
                sprintf(path0,"test/train_190102/%.2f_%.2f_%.2f", alpha, bin_size,max_dist);



                string command ="mkdir -p "+ string(path0); 
                system(command.c_str());

                string path1 = string(path0) + "/count.dat";
                string path2 = string(path0) + "/energy.dat";

                if (misc::file_exists(path1)){
                    continue;
                }

                FILE * count_out =fopen(path1.c_str(),"w");
                FILE * etable_out = fopen(path2.c_str(),"w");

                dfire_calculator * dc = new dfire_calculator();

                dc->train_dfire(pdbl, alpha, bin_size, max_dist);

                dc->print_count(count_out);
                dc->print_edfire(etable_out);

                fclose(count_out);
                fclose(etable_out);

                delete dc;

                cout << path0 << " finished!" << endl;

            }
        }
    }

    // vector<double> vmax_dist ={17.0,18.0,19.0,20.0};    

}


void train_dis_190117(){
    string pdbl = "dataset/PDB/RNA_hires_xp/RNA_hires_xp_pdbs.trainpdbs";

    // vector<double> valpha = {1.30, 1.50, 1.61, 1.70, 1.80, 2.0};
    vector<double> valpha = {1.3, 1.4, 1.5, 1.61, 1.7, 1.8};
    // vector<double> vbin_size = {0.3, 0.5, 1.0};
    // vector<double> vbin_size = {0.5};
    vector<double> vbin_size = {0.3, 0.5, 0.7};
    vector<double> vmax_dist ={15.0, 17.0, 19.0};
    // vector<double> vmax_dist ={15.0, 18.0, 20};
    // vector<double> vmax_dist ={15.0};


    for (auto& alpha: valpha){
        for (auto& bin_size: vbin_size){
            for (auto& max_dist: vmax_dist){

                char path0[200];
                sprintf(path0,"test/train_190117/%.2f_%.2f_%.2f", alpha, bin_size,max_dist);



                string command ="mkdir -p "+ string(path0); 
                system(command.c_str());

                string path1 = string(path0) + "/count.dat";
                string path2 = string(path0) + "/energy.dat";

                if (misc::file_exists(path1)){
                    continue;
                }

                FILE * count_out =fopen(path1.c_str(),"w");
                FILE * etable_out = fopen(path2.c_str(),"w");

                dfire_calculator * dc = new dfire_calculator();

                dc->train_dfire(pdbl, alpha, bin_size, max_dist);

                dc->print_count(count_out);
                dc->print_edfire(etable_out);

                fclose(count_out);
                fclose(etable_out);

                delete dc;

                cout << path0 << " finished!" << endl;

            }
        }
    }

    // vector<double> vmax_dist ={17.0,18.0,19.0,20.0};    

}



void train_dis(){
    string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.l";
    // string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.tail5";
    // string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.tail100";

    vector<double> valpha = {1.61};
    // vector<double> vbin_size = {0.3,0.5,0.7,1.0};
    vector<double> vbin_size = {1.0};
    // vector<double> vmax_dist ={15.0,16.0,17.0,18.0,19.0,20.0};
    vector<double> vmax_dist ={17.0,18.0,19.0,20.0};
    for (auto& alpha: valpha){
        for (auto& bin_size: vbin_size){
            for (auto& max_dist: vmax_dist){

                char path0[200];
                sprintf(path0,"/home/tc/GIT/dfire_rna2/test/train_params/%.2f_%.2f_%.2f", alpha, bin_size,max_dist);

                string command ="mkdir -p "+ string(path0); 
                system(command.c_str());



                // sprintf(path2,"/home/tc/GIT/dfire_rna2/test/train_params/%.2f_%.2f_%.2f/energy.dat", alpha, bin_size,max_dist);
                // system("mkdir -p /home/tc/GIT/dfire_rna2/test/train_params/%.2f_%.2f_%.2f/");

                string path1 = string(path0) + "/count.dat";
                string path2 = string(path0) + "/energy.dat";

                // cout << path1 << endl;

                // cout << path2 << endl;

                FILE * count_out =fopen(path1.c_str(),"w");
                FILE * etable_out = fopen(path2.c_str(),"w");

                dfire_calculator * dc = new dfire_calculator();

                dc->train_dfire(pdbl, alpha, bin_size, max_dist);

                dc->print_count(count_out);
                dc->print_edfire(etable_out);

                fclose(count_out);
                fclose(etable_out);


                delete dc;

            }
        }
    }



    // FILE * count_out =fopen("/home/tc/GIT/dfire_rna2/test/train_params/count.dat","w");
    // FILE * etable_out = fopen("/home/tc/GIT/dfire_rna2/test/train_params/energy.dat","w");

    // // dfire_dihedral * dc = new dfire_dihedral();
    // dfire_calculator * dc = new dfire_calculator();




    // char sout[200];

    // // sprintf(sout, "%.2f_%.2f_%.2f",alpha,bin_size, max_dist);

    // // string myparam = "" + to_string(bin_size) +"_"+ to_string(max_dist);

    // cout << sout << endl;
    


}


// void train_dih(){

//     string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.l";
//     // string etable_out = "/home/tc/GIT/dfire_rna2/test/trainxp/energy.dih.dat";
//     // string count_out = "/home/tc/GIT/dfire_rna2/test/trainxp/count.dih.dat";

//     FILE * count_out =fopen("/home/tc/GIT/dfire_rna2/test/trainxp/count.dih.dat","w");
//     FILE * etable_out = fopen("/home/tc/GIT/dfire_rna2/test/trainxp/energy.dih.dat","w");

//     dfire_dihedral * dc = new dfire_dihedral();

//     dc->train_dfire_dih(pdbl);

//     // dc->print_edfire_dih(stdout);
//     // dc->print_count_dih(stdout);
//     dc->print_edfire_dih(etable_out);
//     dc->print_count_dih(count_out);

//     fclose(count_out);
//     fclose(etable_out);
// }

// void train_dipolar(){

//     // string pdbl = "/home/tc/GIT/dfire_rna2/test/pdb.l";

//     string pdbl = "/home/tc/GIT/dfire_rna2/dataset/PDB/RNA_hires_xp/rmdecoy.l";

//     string count_pn = "/home/tc/GIT/dfire_rna2/test/trainxp/count.dipo.pn.dat";
//     string count_pp = "/home/tc/GIT/dfire_rna2/test/trainxp/count.dipo.pp.dat";
//     string etable_pn = "/home/tc/GIT/dfire_rna2/test/trainxp/energy.dipo.pn.dat";
//     string etable_pp = "/home/tc/GIT/dfire_rna2/test/trainxp/energy.dipo.pp.dat";

//     dfire_dipolar * dc = new dfire_dipolar();



//     dc->train_dfire_dipo(pdbl, 500);

//     // dc->print_count_dipo(stdout, stdout);

//     dc->print_count_dipo(count_pn, count_pp);
//     dc->print_edfire_dipo(etable_pn, etable_pp);
// }






int main(int argc, char const *argv[])
{
    // train_dis_190102();
    train_dis_190117();
    // train_dis();
    // train_3drna();
    // train_dipolar();
    /* code */
    return 0;
}
