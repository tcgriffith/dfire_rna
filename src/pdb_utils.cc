#include "pdb_utils.h"

namespace pdb_utils{

    double hbond_cutoff=3.2; // hydrogen bond cutoff

// ==============================================
// smcra/RNA
// special RNA NTs
// ==============================================

    std::unordered_map<std::string,std::string> rna_specials({
        {"A", "A"},
        {"U", "U"},
        {"G", "G"},
        {"C", "C"},
        {"1MA", "A" },
        {"1MG", "G" },
        {"2MG", "G" },
        {"5BU", "U" },
        {"5MC", "C" },
        {"5MU", "U" },
        {"6IA", "A" },
        {"7MG", "G" },
        {"A23", "A" },
        {"CCC", "C" },
        {"FHU", "U" },
        {"FMU", "U" },
        {"H2U", "U" },
        {"M2G", "G" },
        {"OMC", "C" },
        {"OMG", "G" },
        {"PPU", "U" },
        {"PSU", "U" },
        {"UD5", "U" },
        {"YG", "G" },  
    });


// ==============================================
// smcra/RNA
// Deal with old/new PDB conventions
// ==============================================

    std::vector<std::string> rna_sideatom({
        "A_C2",
        "A_C4",
        "A_C5",
        "A_C6",
        "A_C8",
        "A_N1",
        "A_N3",
        "A_N6",
        "A_N7",
        "A_N9",
        "C_C2",
        "C_C4",
        "C_C5",
        "C_C6",
        "C_N1",
        "C_N3",
        "C_N4",
        "C_O2",
        "G_C2",
        "G_C4",
        "G_C5",
        "G_C6",
        "G_C8",
        "G_N1",
        "G_N2",
        "G_N3",
        "G_N7",
        "G_N9",
        "G_O6",
        "U_C2",
        "U_C4",
        "U_C5",
        "U_C6",
        "U_N1",
        "U_N3",
        "U_O2",
        "U_O4"
    });

    bool is_rna_sideatom(const std::string& atype_s){
        bool flag = false;

        if (std::find(rna_sideatom.begin(), rna_sideatom.end(), atype_s) != rna_sideatom.end()){
            flag = true;
        }

        return flag;
    } 


    std::unordered_map<std::string,std::string> old_2_new({
        {"C1*", "C1'"},
        {"C2*", "C2'"},
        {"C3*", "C3'"},
        {"C4*", "C4'"},
        {"C5*", "C5'"},
        {"O2*", "O2'"},
        {"O3*", "O3'"},
        {"O4*", "O4'"},
        {"O5*", "O5'"},
        {"O1P", "OP1"},
        {"O2P", "OP2"},
        {"rC","C"},
        {"rU","U"},
        {"rG","G"},
        {"rA","A"},

        // add GTP -> G
        {"GTP", "G"},
        {"GDP", "G"},
    });




// ==============================================
// neighbor atom maps in RNA
// ==============================================

    std::unordered_map<std::string, std::vector<std::string>> rna_atom_neighbors({
        {"A_N1",{"C2","C6"}},
        {"A_N3",{"C2","C4"}},
        {"A_N6",{"C6"}},
        {"A_N7",{"C5","C8"}},
        {"A_OP1",{"P"}},
        {"A_OP2",{"P"}},
        {"C_N3",{"C2","C4"}},
        {"C_N4",{"C4"}},
        {"C_O2",{"C2"}},
        {"C_OP1",{"P"}},
        {"C_OP2",{"P"}},
        {"G_N1",{"C2","C6"}},
        {"G_N2",{"C2"}},
        {"G_N3",{"C2","C4"}},
        {"G_N7",{"C5","C8"}},
        {"G_O6",{"C6"}},
        {"G_OP1",{"P"}},
        {"G_OP2",{"P"}},
        {"U_N3",{"C2","C4"}},
        {"U_O2",{"C2"}},
        {"U_O4",{"C4"}},
        {"U_OP1",{"P"}},
        {"U_OP2",{"P"}},
        {"A_O2'",{"C2'"}},
        {"A_O3'",{"C3'","P_n"}},
        {"A_O4'",{"C1'","C4'"}},
        {"A_O5'",{"P","C5'"}},
        {"C_O2'",{"C2'"}},
        {"C_O3'",{"C3'","P_n"}},
        {"C_O4'",{"C1'","C4'"}},
        {"C_O5'",{"P","C5'"}},
        {"G_O2'",{"C2'"}},
        {"G_O3'",{"C3'","P_n"}},
        {"G_O4'",{"C1'","C4'"}},
        {"G_O5'",{"P","C5'"}},
        {"U_O2'",{"C2'"}},
        {"U_O3'",{"C3'","P_n"}},
        {"U_O4'",{"C1'","C4'"}},
        {"U_O5'",{"P","C5'"}},
    });


// ==============================================
// Protein Triple code to single char code
// ==============================================


    std::unordered_map<std::string, std::string> prot_tri2sin({
        {"ALA","A"},
        {"CYS","C"},
        {"ASP","D"},
        {"GLU","E"},
        {"PHE","F"},
        {"GLY","G"},
        {"HIS","H"},
        {"ILE","I"},
        {"LYS","K"},
        {"LEU","L"},
        {"MET","M"},
        {"ASN","N"},
        {"PRO","P"},
        {"GLN","Q"},
        {"ARG","R"},
        {"SER","S"},
        {"THR","T"},
        {"VAL","V"},
        {"TRP","W"},
        {"TYR","Y"},
        {"SEC","U"},
        {"PYL","O"},
        {"ASX","B"},
        {"GLX","Z"},
        {"XLE","J"},
        {"XAA","X"},
    });


//  ==============================================
//  RNA polar atom type encoding
//  RNA {atom name,int}  <-> {int,atom name}
//             polartype <-> typepolar
//  ==============================================

    std::unordered_map<std::string, int> polartype({
        {"A_N1",0},
        {"A_N3",1},
        {"A_N6",2},
        {"A_N7",3},
        {"A_OP1",4},
        {"A_OP2",5},
        {"C_N3",6},
        {"C_N4",7},
        {"C_O2",8},
        {"C_OP1",9},
        {"C_OP2",10},
        {"G_N1",11},
        {"G_N2",12},
        {"G_N3",13},
        {"G_N7",14},
        {"G_O6",15},
        {"G_OP1",16},
        {"G_OP2",17},
        {"U_N3",18},
        {"U_O2",19},
        {"U_O4",20},
        {"U_OP1",21},
        {"U_OP2",22},
        {"A_O2'",23},
        {"A_O3'",24},
        {"A_O4'",25},
        {"A_O5'",26},
        {"C_O2'",27},
        {"C_O3'",28},
        {"C_O4'",29},
        {"C_O5'",30},
        {"G_O2'",31},
        {"G_O3'",32},
        {"G_O4'",33},
        {"G_O5'",34},
        {"U_O2'",35},
        {"U_O3'",36},
        {"U_O4'",37},
        {"U_O5'",38},
    });
    std::unordered_map<int, std::string> typepolar({
        {0,"A_N1"},
        {1,"A_N3"},
        {2,"A_N6"},
        {3,"A_N7"},
        {4,"A_OP1"},
        {5,"A_OP2"},
        {6,"C_N3"},
        {7,"C_N4"},
        {8,"C_O2"},
        {9,"C_OP1"},
        {10,"C_OP2"},
        {11,"G_N1"},
        {12,"G_N2"},
        {13,"G_N3"},
        {14,"G_N7"},
        {15,"G_O6"},
        {16,"G_OP1"},
        {17,"G_OP2"},
        {18,"U_N3"},
        {19,"U_O2"},
        {20,"U_O4"},
        {21,"U_OP1"},
        {22,"U_OP2"},
        {23,"A_O2'"},
        {24,"A_O3'"},
        {25,"A_O4'"},
        {26,"A_O5'"},
        {27,"C_O2'"},
        {28,"C_O3'"},
        {29,"C_O4'"},
        {30,"C_O5'"},
        {31,"G_O2'"},
        {32,"G_O3'"},
        {33,"G_O4'"},
        {34,"G_O5'"},
        {35,"U_O2'"},
        {36,"U_O3'"},
        {37,"U_O4'"},
        {38,"U_O5'"},
    });

    


//  ==============================================
//  RNA atom type encoding
//  RNA {atom name,int} <-> {int,atom name}
//             atomtype <-> typeatom
//  ==============================================
    std::unordered_map<std::string,int> atomtype({
        {"A_C1'", 0},
        {"A_C2",1},
        {"A_C2'",2},
        {"A_C3'",3},
        {"A_C4",4},
        {"A_C4'",5},
        {"A_C5",6},
        {"A_C5'",7},
        {"A_C6",8},
        {"A_C8",9},
        {"A_N1",10},
        {"A_N3",11},
        {"A_N6",12},
        {"A_N7",13},
        {"A_N9",14},
        {"A_O2'",15},
        {"A_O3'",16},
        {"A_O4'",17},
        {"A_O5'",18},
        {"A_OP1",19},
        {"A_OP2",20},
        {"A_P",21},
        {"C_C1'",22},
        {"C_C2",23},
        {"C_C2'",24},
        {"C_C3'",25},
        {"C_C4",26},
        {"C_C4'",27},
        {"C_C5",28},
        {"C_C5'",29},
        {"C_C6",30},
        {"C_N1",31},
        {"C_N3",32},
        {"C_N4",33},
        {"C_O2",34},
        {"C_O2'",35},
        {"C_O3'",36},
        {"C_O4'",37},
        {"C_O5'",38},
        {"C_OP1",39},
        {"C_OP2",40},
        {"C_P",41},
        {"G_C1'",42},
        {"G_C2",43},
        {"G_C2'",44},
        {"G_C3'",45},
        {"G_C4",46},
        {"G_C4'",47},
        {"G_C5",48},
        {"G_C5'",49},
        {"G_C6",50},
        {"G_C8",51},
        {"G_N1",52},
        {"G_N2",53},
        {"G_N3",54},
        {"G_N7",55},
        {"G_N9",56},
        {"G_O2'",57},
        {"G_O3'",58},
        {"G_O4'",59},
        {"G_O5'",60},
        {"G_O6",61},
        {"G_OP1",62},
        {"G_OP2",63},
        {"G_P",64},
        {"U_C1'",65},
        {"U_C2",66},
        {"U_C2'",67},
        {"U_C3'",68},
        {"U_C4",69},
        {"U_C4'",70},
        {"U_C5",71},
        {"U_C5'",72},
        {"U_C6",73},
        {"U_N1",74},
        {"U_N3",75},
        {"U_O2",76},
        {"U_O2'",77},
        {"U_O3'",78},
        {"U_O4",79},
        {"U_O4'",80},
        {"U_O5'",81},
        {"U_OP1",82},
        {"U_OP2",83},
        {"U_P",84},
    });

    std::unordered_map<int, std::string> typeatom({
        {0,"A_C1'"},
        {1,"A_C2"},
        {2,"A_C2'"},
        {3,"A_C3'"},
        {4,"A_C4"},
        {5,"A_C4'"},
        {6,"A_C5"},
        {7,"A_C5'"},
        {8,"A_C6"},
        {9,"A_C8"},
        {10,"A_N1"},
        {11,"A_N3"},
        {12,"A_N6"},
        {13,"A_N7"},
        {14,"A_N9"},
        {15,"A_O2'"},
        {16,"A_O3'"},
        {17,"A_O4'"},
        {18,"A_O5'"},
        {19,"A_OP1"},
        {20,"A_OP2"},
        {21,"A_P"},
        {22,"C_C1'"},
        {23,"C_C2"},
        {24,"C_C2'"},
        {25,"C_C3'"},
        {26,"C_C4"},
        {27,"C_C4'"},
        {28,"C_C5"},
        {29,"C_C5'"},
        {30,"C_C6"},
        {31,"C_N1"},
        {32,"C_N3"},
        {33,"C_N4"},
        {34,"C_O2"},
        {35,"C_O2'"},
        {36,"C_O3'"},
        {37,"C_O4'"},
        {38,"C_O5'"},
        {39,"C_OP1"},
        {40,"C_OP2"},
        {41,"C_P"},
        {42,"G_C1'"},
        {43,"G_C2"},
        {44,"G_C2'"},
        {45,"G_C3'"},
        {46,"G_C4"},
        {47,"G_C4'"},
        {48,"G_C5"},
        {49,"G_C5'"},
        {50,"G_C6"},
        {51,"G_C8"},
        {52,"G_N1"},
        {53,"G_N2"},
        {54,"G_N3"},
        {55,"G_N7"},
        {56,"G_N9"},
        {57,"G_O2'"},
        {58,"G_O3'"},
        {59,"G_O4'"},
        {60,"G_O5'"},
        {61,"G_O6"},
        {62,"G_OP1"},
        {63,"G_OP2"},
        {64,"G_P"},
        {65,"U_C1'"},
        {66,"U_C2"},
        {67,"U_C2'"},
        {68,"U_C3'"},
        {69,"U_C4"},
        {70,"U_C4'"},
        {71,"U_C5"},
        {72,"U_C5'"},
        {73,"U_C6"},
        {74,"U_N1"},
        {75,"U_N3"},
        {76,"U_O2"},
        {77,"U_O2'"},
        {78,"U_O3'"},
        {79,"U_O4"},
        {80,"U_O4'"},
        {81,"U_O5'"},
        {82,"U_OP1"},
        {83,"U_OP2"},
        {84,"U_P"},
    });

//  ==============================================
//  functions
//  ==============================================

    bool is_rna_polar_atom(const std::string& atype){
        bool flag = false;

        if (rna_atom_neighbors.find(atype) != pdb_utils::rna_atom_neighbors.end()){
            flag = true;
        }
        return flag;
    }

    int get_polartype_int(std::string atype_s){
        if (polartype.find(atype_s) != polartype.end()){

            return polartype[atype_s];
        }
        else{
            std::cerr << "polartype error: "<<  atype_s << " not found" << std::endl;
            return -1;
        }        
    }

    int get_atomtype_int(std::string atype_s){
        if (atomtype.find(atype_s) != atomtype.end()){

            return atomtype[atype_s];
        }
        else{
            std::cerr << "atomtype error " <<atype_s << " not found" << std::endl;
            return -1;
        }
    }

// pdb format conversion, eg C4* -> C4', see old_2_new
    std::string format_rna(const std::string oldstr){
        std::string newstr = oldstr;
        auto search = old_2_new.find(oldstr);
        if(search != old_2_new.end()){
            newstr = search->second;
        }
        return newstr;
    }

    std::string protein_tri2sin(const std::string &str_tri){
        auto search = prot_tri2sin.find(str_tri);
        if(search != prot_tri2sin.end()){
            return search->second;
        }
        else {
            std::cerr << &str_tri << " not found" << std::endl;
            return "X";
        }

    }

    // 0 Protein
    // 1 RNA
    // 2 other
    int guess_ctype(const std::string& resname){

        int type = -1;

        auto search = prot_tri2sin.find(resname);
        if (search != prot_tri2sin.end()){
            type = 0;
        } 
        else {
            std::unordered_set<std::string> rna_nt = {"A","U","G","C"};

            if (rna_nt.count(resname)){
                type = 1;
            }
            else{
                type = 2;
            }

        }

        return type;
    }

    std::string rna_tri2sin(const std::string &str_tri){
        auto search = rna_specials.find(str_tri);
        if(search != rna_specials.end()){
            return search->second;
        }
        else {
            std::cerr << &str_tri << " not found" << std::endl;
            return "X";
        }        
    }
}
