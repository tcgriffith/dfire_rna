#include "PDB.h"


// ###########################################
// Atom
// ###########################################

double atom::distance(atom &other){
    double dist=0.0;

    double *x1 = this->get_x();
    double *x2 = other.get_x();

    dist= sqrt(pow((x1[0]-x2[0]),2)+
               pow((x1[1]-x2[1]),2)+
               pow((x1[2]-x2[2]),2)
               );

    return dist; 
}

double atom::distance2(atom &other){
    double dist2=0.0;
    double *x1 = this->get_x();
    double *x2 = other.get_x();

    dist2 = (x1[0] - x2[0]) * (x1[0] - x2[0])+
            (x1[1] - x2[1]) * (x1[1] - x2[1])+
            (x1[2] - x2[2]) * (x1[2] - x2[2]);

    // dist2= (pow((x1[0]-x2[0]),2)+
    //            pow((x1[1]-x2[1]),2)+
    //            pow((x1[2]-x2[2]),2)
    //            );

    return dist2;
}

double atom::get_cos_angle_p(atom &other){
    double anglep=-999;

    if (this->is_polar()){
        Xvec p(this->get_x());
        Xvec q(other.get_x());
        Xvec &rref = this->get_pvec();

        Xvec pq = q - p;

        anglep =  rref.angle_cos(pq);
    }

    return anglep;
}

double atom::get_cos_angle_pq(atom &other){
    double anglep=-999;

    if (this->is_polar() && other.is_polar()){

        Xvec &p_ref = this->get_pvec();

        Xvec &q_ref = other.get_pvec();

        anglep =  p_ref.angle_cos(q_ref);
    }

    return anglep;
}

double atom::torsion(atom &aa, atom &ab, atom &ac, atom &ad){
    Xvec x1(aa.get_x());
    Xvec x2(ab.get_x());
    Xvec x3(ac.get_x());
    Xvec x4(ad.get_x());
    double tor = torsion_xvec(x1,x2,x3,x4);
    return tor;
}


// ###########################################
// Residue
// ###########################################


bool residue::rna_check_atoms(){
    bool atomsfull=true;
    // if (atoms_map.find("O3'") == atoms_map.end()) atomsfull=false;
    if (atoms_map.find("P") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("P");}
    if (atoms_map.find("O5'") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("O5'");}
    if (atoms_map.find("C5'") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("C5'");}
    if (atoms_map.find("C4'") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("C4'");}
    if (atoms_map.find("C3'") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("C3'");}
    if (atoms_map.find("O3'") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("O3'");}
    // if (atoms_map.find("P") == atoms_map.end()) atomsfull=false;
    // if (atoms_map.find("O5'") == atoms_map.end()) atomsfull=false;
    if (atoms_map.find("C2'") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("C2'");}
    if (atoms_map.find("C1'") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("C1'");}


    // if (atoms_map.find(an) ==atoms_map.end()) atomsfull=false;

    if (this->get_resname() == "A" || this->get_resname() == "G"){
        if (atoms_map.find("N9") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("N9");}
        if (atoms_map.find("C8") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("C8");}

    }

    if(this->get_resname() == "U" || this->get_resname() == "C"){
        if (atoms_map.find("N1") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("N1");}
        if (atoms_map.find("C6") == atoms_map.end()) {atomsfull=false; this->lack_atoms.push_back("C6");}

    }

    return atomsfull;

}

// check hbond_atoms
bool residue::rna_check_hhatoms(){
    bool atomsfull=true;

    if (this->get_resname() == "A"){
        if (atoms_map.find("N1") == atoms_map.end()) atomsfull=false;
        if (atoms_map.find("N6") == atoms_map.end()) atomsfull=false;
    }
    else if (this->get_resname() == "G"){
        if (atoms_map.find("N2") == atoms_map.end()) atomsfull=false;
        if (atoms_map.find("N1") == atoms_map.end()) atomsfull=false;
        if (atoms_map.find("O6") == atoms_map.end()) atomsfull=false;
        if (atoms_map.find("N7") == atoms_map.end()) atomsfull=false;
    }
    else if (this->get_resname() == "C"){
        if (atoms_map.find("O2") == atoms_map.end()) atomsfull=false;
        if (atoms_map.find("N3") == atoms_map.end()) atomsfull=false;
        if (atoms_map.find("N4") == atoms_map.end()) atomsfull=false;

    }
    else if (this->get_resname() == "U"){
        if (atoms_map.find("N3") == atoms_map.end()) atomsfull=false;
        if (atoms_map.find("O4") == atoms_map.end()) atomsfull=false;
        if (atoms_map.find("O2") == atoms_map.end()) atomsfull=false;
    }



    return atomsfull;
}


// test if residue has atoms to form hydrogen bonds;

std::vector<atom *> residue::rna_get_hhatoms(){
    std::vector<atom *> hhatoms;

    if (this->get_resname() == "A"){
        hhatoms.push_back( this->get_atom("N1"));
        hhatoms.push_back( this->get_atom("N6"));
        hhatoms.push_back( this->get_atom("N7"));
    }
    else if (this->get_resname() == "G"){
        hhatoms.push_back( this->get_atom("N2"));
        hhatoms.push_back( this->get_atom("N1"));
        hhatoms.push_back( this->get_atom("O6"));
        hhatoms.push_back( this->get_atom("N7"));
        hhatoms.push_back( this->get_atom("N3"));
        hhatoms.push_back( this->get_atom("N9"));
    }
    else if (this->get_resname() == "C"){
        hhatoms.push_back( this->get_atom("O2"));
        hhatoms.push_back( this->get_atom("N3"));
        hhatoms.push_back( this->get_atom("N4"));

    }
    else if (this->get_resname() == "U"){
        hhatoms.push_back( this->get_atom("N3"));
        hhatoms.push_back( this->get_atom("O4"));
        hhatoms.push_back( this->get_atom("O2"));
    }

    return hhatoms;

}

Xvec residue::rna_get_planevec(){
    std::vector<atom *> scatoms;
    if (this->get_resname() == "A" || this->get_resname() == "G") {
        scatoms.push_back(get_atom("C1'"));
        scatoms.push_back(get_atom("C8"));
        scatoms.push_back(get_atom("N1"));

    }

    else if(this->get_resname() == "U" || this->get_resname() == "C"){
        scatoms.push_back(get_atom("C1'"));
        scatoms.push_back(get_atom("C6"));
        scatoms.push_back(get_atom("N3"));
    }

    Xvec a = Xvec(scatoms[0]->get_x());
    Xvec b = Xvec(scatoms[1]->get_x());
    Xvec c = Xvec(scatoms[2]->get_x());

    Xvec plane_vec = (c - a).cross_product(b - a);

    return plane_vec;

}



bool residue::is_rna_pair(residue & other_res){

    if (!this->rna_check_hhatoms() || !other_res.rna_check_hhatoms()) return false;

    atom * c4a = this->get_atom("C4'");
    atom * c4b = other_res.get_atom("C4'");

    if (c4a->distance2(*c4b)> 500) return false;

    bool is_pair=false;

    std::vector<atom *> res_a_hhatm = this->rna_get_hhatoms();
    std::vector<atom *> res_b_hhatm = other_res.rna_get_hhatoms();

    // std::vector<atom *> res_a_sc = this->get_scatoms();
    // std::vector<atom *> res_b_sc = other_res.get_scatoms();

    Xvec plane_a = this->rna_get_planevec();

    int n_pair=0;

    for (atom * p_atom_a : res_a_hhatm){
        for ( atom * p_atom_b : res_b_hhatm){
            if (p_atom_a->distance(* p_atom_b) < pdb_utils::hbond_cutoff){

                Xvec vec_hh = Xvec (p_atom_b->get_x()) - Xvec(p_atom_a->get_x());

                double angle_cos = plane_a.angle_cos(vec_hh);
                // fprintf(stdout, "%d %d %s %s %s %s %8.6f\n",this->get_resid(), other_res.get_resid(), this->get_resname().c_str(),p_atom_a->get_name().c_str(),other_res.get_resname().c_str(), p_atom_b->get_name().c_str(), angle_cos);
                //  angle of hydrogen vec and plane vec should sit between 70 ~110
                if (abs(angle_cos) < .342)
                n_pair ++;
            }
        }
    }

    if (n_pair >= 2) is_pair=true;
    return is_pair;
}


//###########################################
// Chain
//###########################################





std::string chain::get_seq(){

    std::string seq="";

    for (int i = 0; i < this->get_residue_num(); ++i){
        seq=seq+this->residues[i].get_ressin();
    }

    // if (this->chaintype == "PROTEIN"){
    //     for (int i = 0; i < this->get_residue_num(); ++i){
    //         seq=seq+this->residues[i].get_ressin();
    //     }
    // }
    // else {
    //     for (int i = 0; i < this->get_residue_num(); ++i){
    //         seq=seq+this->residues[i].get_resname();
    //     }
    // }

    return seq;

}




// ###########################################
// Structure
// ###########################################


bool structure::readpdb(std::string pdbname){

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

    chain *a_chain;
    residue *a_res;
    atom *a_atom;

    // int nchain,nres,natom;

    nchain=nres=natom=0;



    while(getline(fs,line)){
        if(line.substr(0,3) == "END") break;
        if(line.substr(0,3) == "TER") continue;
        if(line.substr(0,4) != "ATOM") continue;

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
                this->add_chain(*(new chain()));
                a_chain =&(this->chains.back());
                a_chain->set_chainid(cid);
                // guess chain type 
                int type = pdb_utils::guess_ctype(rn);

                if (type == 0) a_chain->set_chaintype("PROTEIN");
                else if (type == 1) a_chain->set_chaintype("NT");
                lastcid = cid;
                nchain++;
            }

            a_chain->add_residue(*(new residue()));
            a_res = &(a_chain->residues.back());
            a_res->set_name(rn);
            a_res->set_resid(rid);
        // make residue remember chainid
            a_res->set_chainid(cid);
        // set ressin for Protein and RNA
            if (a_chain->get_chaintype() == "PROTEIN") {
                string res_sin = pdb_utils::protein_tri2sin(rn);
                a_res->set_ressin(res_sin);
            }
            if (a_chain->get_chaintype() == "NT") {
                string res_sin = rn;
                if (rn.length() > 1){
                    res_sin = rn.back();
                    cerr << "# warning: unknown RNA NT type: " << rn << " set as:" << res_sin << endl;
                }

                a_res->set_ressin(res_sin);
            }

            lastrid = rid;
            nres++;
        }

        a_atom = new atom();

        
        // a_atom = &(a_res->atoms.back());
        scrd=line.substr(30,24);
        a_atom->set_x(scrd);
        a_atom->set_name(an);
        a_atom->set_type(rn+"_"+an);

        a_res->addatom( *a_atom);

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

bool structure::rna_init_pairs(){

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

bool structure::rna_init_polar(){

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
                if (pdb_utils::is_rna_polar_atom(atype)) {
                    a_atom.set_polar(true);
 
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
                            cout << i << " error, missing atom: " << nb_atom_name << " " << a_res->get_resid() << endl;
                            return 1;
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

std::string structure::get_mapstr(){

    std::string outstr = "";
    char tmp[2000];

    sprintf(tmp,"LEN %d\n", this->get_nres());

    outstr += tmp;


    for (auto & apair: this->pairs){
        sprintf(tmp,"CON\t%s\t%s\t1\n",apair.first->get_resid().c_str(), apair.second->get_resid().c_str());

        outstr += tmp;
    }

    return outstr;
}
