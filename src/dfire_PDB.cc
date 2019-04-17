#include "dfire_PDB.h"


// dfire specific
// replaced:
// chain -> rnachain
// residue -> rnabase
// structure -> rna




// test if residue has atoms to form hydrogen bonds;
 


void rnabase::calc_torsion(rnabase * r_up, rnabase * r_down){

    atom * main_chain[9];
    atom * side_chain[4];

    
// r_up ,this and rdown should be adjacent on sequence

   // if (r_up->get_resid() == this->get_resid() -1
   //      && r_down->get_resid() == this->get_resid() +1
   //      && r_up->check_atoms() && r_down->check_atoms() && this->check_atoms()){
    if (r_up->get_residsd() == this->get_residsd() -1
        && r_down->get_residsd() == this->get_residsd() +1
        && r_up->rna_check_atoms() && r_down->rna_check_atoms() && this->rna_check_atoms()){
        main_chain[0]=r_up->get_atom("O3'");
        main_chain[1]=get_atom("P");
        main_chain[2]=get_atom("O5'");
        main_chain[3]=get_atom("C5'");
        main_chain[4]=get_atom("C4'");
        main_chain[5]=get_atom("C3'");
        main_chain[6]=get_atom("O3'"); 
        main_chain[7]=r_down->get_atom("P");
        main_chain[8]=r_down->get_atom("O5'");

        side_chain[0]=get_atom("C2'");
        side_chain[1]=get_atom("C1'");

        if (this->get_resname() == "A" || this->get_resname() == "G"){
            side_chain[2]=get_atom("N9");
            side_chain[3]=get_atom("C8");
        }
        else if(this->get_resname() == "U" || this->get_resname() == "C"){
            side_chain[2]=get_atom("N1");
            side_chain[3]=get_atom("C6");
        }


    // alpha
        this->dihedrals[0]= atom::torsion(* main_chain[0], * main_chain[1], * main_chain[2], * main_chain[3]);
    // beta
        this->dihedrals[1]= atom::torsion(* main_chain[1], * main_chain[2], * main_chain[3], * main_chain[4]);
    // gamma
        this->dihedrals[2]= atom::torsion(* main_chain[2], * main_chain[3], * main_chain[4], * main_chain[5]);
    // delta
        this->dihedrals[3]= atom::torsion(* main_chain[3], * main_chain[4], * main_chain[5], * main_chain[6]);
    // epsilon
        this->dihedrals[4]= atom::torsion(* main_chain[4], * main_chain[5], * main_chain[6], * main_chain[7]);
    // zeta
        this->dihedrals[5]= atom::torsion(* main_chain[5], * main_chain[6], * main_chain[7], * main_chain[8]);
    // chi
        this->dihedrals[6]= atom::torsion(* side_chain[0], * side_chain[1], * side_chain[2], * side_chain[3]);



        this->_has_dih=true;

   }
   else {
        this->_has_dih=false;
   }
}

bool rnachain::init_dihedrals(){

    // printf("residues to init: %d",this->get_residue_num()); 
    for (int i = 1; i < this->residues.size() - 1; ++i){
        this->residues[i].calc_torsion(&residues[i-1], &residues[i+1]);

        // // if (this->residues[i].has_dih())
        // {
        //     printf("residue %d initialized\n", i );
        // }
    }

    return 0;
}


bool rna::readpdb(std::string pdbname){

// string filepath="/home/tc/GIT/dfire_rna2/test/5o62.pdb";

// rid should be string, see 3g8t.pdb, line 1779

    ifstream fs(pdbname);

    string line;
   // residue atom chain name;
    string rn, an;
    string lastrid="RAGNAROK!";
    string lastcid="RAGNAROK!";
    string lastalt="RAGNAROK!";
    string lastalt_res = "KARAGE!";
    string rid;
    string cid;
    string scrd;

    string alt;
    bool has_alt = false;

    rnachain *a_chain;
    rnabase *a_res;
    atom *a_atom;

    // int nchain,nres,natom;

    nchain=nres=natom=0;



    while(getline(fs,line)){
        if(line.substr(0,3) == "END") break;
        if(line.substr(0,3) == "TER") continue;
        if(line.substr(0,4) != "ATOM") continue;
    // ignore H atom
        if(line.find("H") != std::string::npos) continue; 

        rn = line.substr(17,3);
        rn = misc::trim(rn);
        rn = pdb_utils::format_rna(rn);
        an = line.substr(13,3);
        an = pdb_utils::format_rna(an);

        alt = line.substr(16,1);
        an = misc::trim(an); 
        cid=line[21];
        rid = line.substr(22,5);
        rid = misc::trim(rid);

        if (alt != " "){
            if (lastalt == "RAGNAROK!") lastalt = alt;
            else if (alt != lastalt) {
                has_alt = true;
                continue;
            }

        }



        if (rid != lastrid || cid != lastcid){

            if (cid != lastcid){
                // cout << cid <<endl;
                // this->add_chain(*(new rnachain()));
                this->chains.push_back(*(new rnachain()));
                a_chain =&(this->chains.back());
                a_chain->set_chainid(cid);
                // guess chain type by res name length
                if (rn.length() ==3) a_chain->set_chaintype("PROTEIN");
                else a_chain->set_chaintype("NT");
                lastcid = cid;

                nchain++;
            }
                
            // a_chain->add_residue(*(new residue()));
            a_chain->residues.push_back(*(new rnabase()));
            a_res = &(a_chain->residues.back());
            a_res->set_name(rn);
            a_res->set_resid(rid);
        // make residue remember chainid
            a_res->set_chainid(cid);
        // set ressin for chain type = PROTEIN
            if (a_chain->get_chaintype() == "PROTEIN") {
                string res_sin = pdb_utils::protein_tri2sin(rn);

                a_res->set_ressin(res_sin);

            }

            lastrid = rid;
            nres++;
        }

        a_res->addatom(*(new atom()));
        a_atom = &(a_res->atoms.back());
        scrd=line.substr(30,24);
        a_atom->set_x(scrd);
        a_atom->set_name(an);
        a_atom->set_type(rn+"_"+an);
        if (pdb_utils::is_rna_polar_atom(a_atom->get_type())) a_atom->set_polar(true);

        natom++;
    }

    fs.close();

    int ridsd =0;

    for(auto & achain: this->chains) {
        ridsd =0;
        for ( auto & aresidue : achain.residues){
            aresidue.set_residsd(ridsd++);
            for (auto & aatom : aresidue.atoms ){
                aresidue.addatom_tomap( &aatom, aatom.get_name());
            }
        }
    }
    if (has_alt) {
        cerr << "# Warning: alternative conformation detected. First alt conformation used; File: " << pdbname  << endl;

    }

    return true;
}


// override structure 
bool rna::rna_init_pairs(){

    if (!this->_polar_inited) rna_init_polar();



    for (auto & achain : this->chains){
        for (auto & bchain: this->chains){
            for (auto & ares : achain.residues){
                for (auto & bres : bchain.residues){


                    if (achain.get_chainid() > bchain.get_chainid()) continue;
                    if (achain.get_chainid() == bchain.get_chainid()){
                        if (ares.get_residsd() >= bres.get_residsd()) continue;
                    }

                    if (ares.is_rna_pair(bres)){

                        this->pairs.push_back(pair<residue *, residue *>(&ares,&bres));

                    }

                }
            }
        }
    }

    this->_pair_inited = true;

    return 0;
}
 
bool rna::rna_init_polar(){

    // unordered_set<char> polar_set = {'O','N','P'};
    for(auto & achain: this->chains){
        for (int i=0; i < achain.residues.size(); i++){


            residue * a_res = &achain.residues.at(i);
            residue * next_res;
    // residues should pass atoms check
            if (!a_res->rna_check_atoms()) continue;

            if (i + 1 < achain.residues.size()) next_res = &achain.residues.at(i+1);
            else next_res =nullptr;


            for (auto & a_atom : a_res->atoms){
                string atype = a_atom.get_type();
                if (a_atom.is_polar()){
                // if (pdb_utils::is_rna_polar_atom(atype)) {
                    // a_atom.set_polar(true);

                    Xvec r_ref;

                    Xvec a_atom_vec = Xvec(a_atom.get_x());

                    for (auto & nb_atom_name : pdb_utils::rna_atom_neighbors.at(atype)){

                        atom * nb_atom;
                        if (nb_atom_name == "P_n"){
                            if (next_res != nullptr && next_res->rna_check_atoms()){
                                nb_atom = next_res->get_atom("P");
                            }
                        }
                        else {
                            nb_atom = a_res->get_atom(nb_atom_name);
                        }

                        if (nb_atom == nullptr){
                            cerr << i << " error, missing atom: " << nb_atom_name << " at residue: " << a_res->get_resid() << endl;
                            return(1);
                        }

                        Xvec nb_atom_vec = Xvec(nb_atom->get_x());
                        r_ref = r_ref + a_atom_vec - nb_atom_vec;
                    }

                    a_atom.set_pvec(r_ref);
                }
            }
        }
    }
    this->_polar_inited = true;
    return 0;
}