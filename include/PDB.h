// no boost this time
// class: atom, residue, chain, structure
// framework for residues
#ifndef _PDB
#define _PDB
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include "xvec.h"
#include "misc.h"
#include "pdb_utils.h"


class atom;
class residue;
class chain;
class structure;


class atom
{
protected:
    double x[3]; // coordinates
    std::string atomname; // atom name from pdb
    std::string atomtype; // atom type from pdb
    int atomid; // atom id from pdb
//  info
    bool _is_polar = false; // based on atom type, O N P -> polar
    std::vector<atom * > nb_atoms;
    Xvec _pvec;

public:
    // basic
    atom(){ atomtype=""; atomname="";atomid=-1;for (int m=0;m <3;m++) x[m]=-1;};
    ~atom(){};
    void set_name(std::string aname){atomname=aname;}
    void set_type(std::string atype){atomtype = atype;}
    void set_id(int aid){atomid = aid;}
    void set_x(const double *xt) {for (int m=0;m <3;m++) x[m]=xt[m];}
    bool set_x(string scrd){
        istringstream iss(scrd);
        for(int i=0; i<3; i++){
            if( (iss>>x[i]).fail() )  return false;
        }
        return true;
    }
    double *get_x(){return x;};
    std::string get_name(){return atomname;};
    std::string get_type(){return atomtype;};
    int get_id(){return atomid;}
    double distance(atom &other);
    double distance2(atom &other);
    static double torsion(atom &aa, atom &ab, atom &ac, atom &ad);


    // polar atoms

    void set_polar(bool polar) {_is_polar = polar;}
    bool is_polar() {return _is_polar;}
    bool set_pvec(Xvec &pv){_pvec = pv; return 0;}
    Xvec & get_pvec(){return _pvec;}
    double get_cos_angle_p(atom &other);
    double get_cos_angle_p(atom * otherp) {return get_cos_angle_p(*otherp);}
    double get_cos_angle_pq(atom &other);
    double get_cos_angle_pq(atom * otherp){return get_cos_angle_pq(*otherp);}

};

class residue
{
protected:
    std::string chainid;
    std::string resid;
    int residsd;
    std::string resname;
    std::string restri;
    std::string ressin;
    int resint;
    std::string ss;

    bool _miss_atom=false;

public:
    std::vector <atom> atoms={};
    std::unordered_map<string,atom *> atoms_map={};
    // lack_atoms initialized by check_atoms
    std::vector<string> lack_atoms = {};
    residue(){resid="";resname="";restri="";ressin='-';resint=-1; residsd =-1;};
    ~residue(){};

    void set_name(std::string rn) {resname=rn;};
    std::string get_resname(){return resname;};

    void set_resid(std::string rid){resid=rid;}
    std::string get_resid(){return resid;};

    void set_residsd(int rid){residsd=rid;}
    int get_residsd(){return residsd;};

    void set_resint(int resin){resint=resin;}
    int get_resint() {return resint;};

    void set_ressin(string res_sin){ressin = res_sin;};
    std::string get_ressin(){return ressin;};

    void set_chainid(std::string cid){chainid = cid;};
    std::string get_chainid(){return chainid;};

    void set_ss(std::string a_ss){ss=a_ss;};
    std::string get_ss() {return ss;};

    // virtual bool check_atoms(std::string name);

    // hydrogen atoms

    bool addatom(atom &a_atom){atoms.push_back(a_atom); return 0;};
    bool addatom_tomap(atom *a_atom, string atomname){atoms_map.insert({atomname,a_atom}); return 0;};
    atom * get_atom(string an){
        if (atoms_map.find(an) !=atoms_map.end()) return atoms_map.at(an);
        else {
            cerr << "atom not found: " << an << endl;
            return NULL;
        }

    };
    virtual int atom_num(){return atoms.size();}



    // RNA specific
    bool _paired = false;
    // move to
    // void set_paired(bool paired){_paired = paired;};

    virtual bool is_rna_pair(residue &);
    virtual bool rna_check_atoms();
    virtual bool rna_check_hhatoms();

    std::vector<atom *> rna_get_hhatoms();
    virtual Xvec rna_get_planevec();

};


class chain
{
protected:
    // PROTEIN, DNA/RNA
    string chaintype="";
    std::string chainid;

public:
    std::vector<residue> residues;
    chain(){this->chainid="RAGNAROK";};
    ~chain(){};

    void set_chainid(std::string cid){chainid = cid;}
    void set_chaintype(std::string ctype){chaintype=ctype;}
    virtual int get_residue_num(){return residues.size();}
    bool add_residue(residue a_res){residues.push_back(a_res); return 0;};

    std::string get_seq();
    std::string get_chainid(){return chainid;}
    std::string get_chaintype(){return chaintype;}
};

class structure
{
protected:
    std::string pdbname;
    int nchain, nres,natom;

// RNA specific
    bool _polar_inited = false;
    bool _pair_inited = false;

public:
    std::vector<chain> chains;
    structure(){};

    structure(std::string pdbfile){
        readpdb(pdbfile);
    };
    ~structure(){};

    void set_pdbname(std::string pname){pdbname=pname;};
    std::string get_pdbname(){return pdbname;}
    int get_nchain(){return nchain;};
    int get_nres(){return nres;};
    int get_natom(){return natom;};

    virtual bool readpdb(std::string);
    bool add_chain(chain &achain){chains.push_back(achain); return 0;};

    // RNA specific
    // use virtual should make dfire_PDB working
    virtual bool rna_init_polar();
    virtual bool rna_init_pairs();
    std::vector<std::pair<residue *, residue *>> pairs;

    virtual std::string get_mapstr();

};


#endif
