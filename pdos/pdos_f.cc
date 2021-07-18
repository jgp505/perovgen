//
// 2005. 8.17 Sung-Hoon Lee : spin-unpolarized case
// 2009.10.26 Sung-Hoon Lee : extension to spin-polarized case
// 2017.07.27 Sung-Hoon Lee : extension to f-orbitals
//
// compile: icpc -O -o pdos pdos.cc
//

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <valarray>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>
#include <complex>

using namespace std;

const double tiny = 1e-6;                         // a tiny number

class PDOS_LIST
{
    public:
        PDOS_LIST(string str) { name = str; }
        string name;
        vector<int> ions;
        vector<string> orbitals;
        vector<string> spins;
};

class PHASE_LIST
{
    public:
        PHASE_LIST(string str) { name = str; }
        string name;
        vector<int> kvec;
        vector<int> band;
        vector<string> orbitals;
        vector<string> spins;
};

class POINT
{
    public:
        double pt[3];
};

class PDOS
{

    public:

        PDOS(int argc, char** argv);
        ~PDOS() { }

        void read_OUTCAR();

        void print_BAND();

        void read_PDOS_LIST();

        void read_PHASE_LIST();

        void read_PROCAR();

        void print_DOS();

        void print_PDOS();

        void print_PHASE();

    protected:

        string file_header;
        map<string, int> orbital_table;
        map<string, int> spin_table;
        double profile_cutoff;
        int folding_period;

        bool PRINT_BAND;
        bool PRINT_DOS;
        bool PRINT_PDOS;
        bool PRINT_PROFILE;
        bool PRINT_PHASE;

        bool open_file(string file_type, FILE*& ifp);

        int nkpts;
        int nbands;
        int nions;
        int ispin;
        double alpha_plus_beta;

        valarray<double> kvec_wgt;
        valarray<POINT> kvec;
        valarray<double> band_en;
        valarray<double> band_occ;
        valarray<int>    band_kvec;

        valarray<double> dos_energy_mesh;
        valarray<double> full_dos;
        valarray<double> occupied_dos;
        valarray<double> accumulated_dos;

        double gaussian_width;
        double e_min;
        double e_max;
        double interval;

        void init_dos();

        vector<PDOS_LIST> pdos_list;

        vector<PHASE_LIST> phase_list;

        void print_usage(char** argv);

        void print_dos(string filename);

        void print_band(string filename, string spin_type);

        void print_band_profile(string filename, valarray<double>& band_wgt);

        valarray<float> orbital_wgt;
        valarray<float> wf_real;
        valarray<float> wf_imag;
        bool phase;

        void calc_dos(valarray<double>& band_wgt, valarray<double>& band_en,
            valarray<double>& band_occ, valarray<int>& band_kvec,
            valarray<double>& kvec_wgt);

        void read_orbital_wgt_from_PROCAR(FILE *ifp, valarray<float>& orbital_wgt,
            valarray<float>& wf_real, valarray<float>& wf_imag,
            int nkpts, int nbands, int nions, int ispin);

        void complex_conversion(double real, double imag, double& rho, double& theta);

        void angle2color(double theta, double& R, double& G, double& B);

        double erfcc(double x);

        template<class T>
            bool read_it_from_str(char* str, char* to_find, char delimiter, T& value);

        template<class T>
            char* read_it_from_file(FILE* ifp, char* to_find, char delimiter, T& value);

        bool read_it_from_str(char* str, char* to_find) {
            int i; return read_it_from_str(str, to_find, 0, i);
        }

        char* read_it_from_file(FILE* ifp, char* to_find) {
            int i; return read_it_from_file(ifp, to_find, 0, i);
        }
} ;

int main(int argc, char **argv)
{
    PDOS job(argc, argv);
    job.read_OUTCAR();
    job.print_BAND();
    job.print_DOS();
    job.read_PDOS_LIST();
    job.read_PHASE_LIST();
    job.read_PROCAR();
    job.print_PDOS();
    job.print_PHASE();
}


PDOS::PDOS(int argc, char **argv)
{
    folding_period = 0x0FFFFFFF;
    gaussian_width = 0.1;
    interval = 0.01;
    profile_cutoff = 0.02;

    PRINT_BAND = false;
    PRINT_DOS = false;
    PRINT_PDOS = false;
    PRINT_PROFILE = false;
    PRINT_PHASE = false;

    vector<bool> arg_read(false, argc);
    for (int i=1; i<argc; i++) {
        string arg = argv[i];
        transform(arg.begin(), arg.end(), arg.begin(), (int(*)(int))tolower);
        if (arg == "band")
            arg_read[i] = PRINT_BAND = true;
        if (arg == "dos")
            arg_read[i] = PRINT_DOS = true;
        if (arg == "pdos")
            arg_read[i] = PRINT_PDOS = true;
        if (arg == "profile")
            arg_read[i] = PRINT_PROFILE = true;
        if (arg == "phase")
            arg_read[i] = PRINT_PHASE = true;
        if (arg.find("period=") == 0) {
            arg.replace(0, 7, "");
            folding_period = atoi(arg.c_str());
            arg_read[i] = true;
        }
        if (arg.find("width=") == 0) {
            arg.replace(0, 6, "");
            gaussian_width = max(atof(arg.c_str()), 1e-5);
            arg_read[i] = true;
        }
        if (arg.find("interval=") == 0) {
            arg.replace(0, 9, "");
            interval = max(atof(arg.c_str()), 1e-3);
            arg_read[i] = true;
        }
        if (arg.find("cutoff=") == 0) {
            arg.replace(0, 7, "");
            profile_cutoff = atof(arg.c_str());
            arg_read[i] = true;
        }
    }

    if (!(PRINT_BAND || PRINT_DOS || PRINT_PDOS || PRINT_PROFILE || PRINT_PHASE))
        print_usage(argv);

    file_header = "";
    if (arg_read[1] == false) {
        file_header = argv[1];
        int i = file_header.find(".");
        int j = file_header.size();
        if (i > 0 && i < j)
            file_header.replace(i, j - i, "");
        arg_read[1] = true;
    }

    for (int i=1; i<argc; i++) {
        if (arg_read[i] == false) {
            fprintf(stderr, "Incorrect %d-th argument: %s\n", i + 1,
                argv[i]);
            print_usage(argv);
        }
    }

    if (PRINT_BAND)
        fprintf(stderr, "BAND ");
    if (PRINT_DOS)
        fprintf(stderr, "DOS ");
    if (PRINT_PDOS)
        fprintf(stderr, "PDOS ");
    if (PRINT_PROFILE)
        fprintf(stderr, "PROFILE ");
    if (PRINT_PHASE)
        fprintf(stderr, "PHASE ");
    fprintf(stderr, "will be calculated.\n");
    if (PRINT_BAND && folding_period != 0x0FFFFFFF)
        fprintf(stderr, "folding_period     = %d\n", folding_period);
    if (PRINT_DOS || PRINT_PDOS || PRINT_PROFILE) {
        fprintf(stderr, "gaussian_width     = %g\n", gaussian_width);
        fprintf(stderr, "mesh space for DOS = %g\n", interval);
    }
    if (PRINT_PROFILE)
        fprintf(stderr, "profile_cutoff     = %g\n", profile_cutoff);
    fprintf(stderr, "\n");

    orbital_table["s"]   = 0;
    orbital_table["py"]  = 1;
    orbital_table["pz"]  = 2;
    orbital_table["px"]  = 3;
    orbital_table["dxy"] = 4;
    orbital_table["dyz"] = 5;
    orbital_table["dz2"] = 6;
    orbital_table["dxz"] = 7;
    orbital_table["dx2"] = 8;
    orbital_table["f-3"] = 9;
    orbital_table["f-2"] = 10;
    orbital_table["f-1"] = 11;
    orbital_table["f0"]  = 12;
    orbital_table["f1"]  = 13;
    orbital_table["f2"]  = 14;
    orbital_table["f3"]  = 15;
    orbital_table["tot"] = 16;

    spin_table["up"] = 0;
    spin_table["down"] = 1;
    spin_table["dn"] = 1;
    spin_table["tot"] = 2;
}


void PDOS::print_usage(char** argv)
{
    fprintf(stderr, "usage: %s %s\n", argv[0], "[OUTCAR or *.voc] "
        "[at least one in band dos pdos profile phase] [options]");
    fprintf(stderr, "option: (folding) period=N\n");
    fprintf(stderr, "        (gaussian) width=X\n");
    fprintf(stderr, "        interval=X\n");
    fprintf(stderr, "        (profile) cutoff=X\n");
    exit(-1);
}


bool PDOS::open_file(string file_type, FILE*& ifp)
{
    string file_name;
    if (file_type == "OUTCAR") {
        file_name = (file_header == "") ? "OUTCAR" : file_header + ".voc";
    }
    else if (file_type == "PROCAR") {
        file_name = (file_header == "") ? "PROCAR" : file_header + ".procar";
    }
    else if (file_type == "LIST") {
        file_name = (file_header == "") ? "LIST" : file_header + ".list";
    }
    else if (file_type == "CHGCAR") {
        file_name = (file_header == "") ? "CHGCAR" : file_header + ".chgcar";
    }
    else if (file_type == "LOCPOT") {
        file_name = (file_header == "") ? "LOCPOT" : file_header + ".locpot";
    }

    fprintf(stderr, "reading %s\n", file_name.c_str());
    if ((ifp = fopen(file_name.c_str(), "r")) == NULL) {
        fprintf(stderr, "fail to open %s; return\n", file_name.c_str());
        return false;
    }
    return true;
}


void PDOS::read_OUTCAR()
{
    FILE* ifp = NULL;
    if (!open_file("OUTCAR", ifp)) {
        fprintf(stderr, "exit\n");
        exit(-1);
    }

    char lineinput[256];
    char *strptr;
    double db;
    bool automatic_k_mesh = false;
    nkpts = 0;

    while (fgets(lineinput, 256, ifp) != NULL) {

        if (read_it_from_str(lineinput, "irreducible k-points:")) {
            automatic_k_mesh = true;
            read_it_from_str(lineinput, "Found", ' ', nkpts);
            kvec_wgt.resize(nkpts);
            read_it_from_file(ifp, "Coordinates               Weight");
            for (int i=0; i<nkpts; i++)
                fscanf(ifp, "%lf%lf%lf%lf", &db, &db, &db, &kvec_wgt[i]);
            kvec_wgt /= kvec_wgt.sum();
        }

        if (! automatic_k_mesh) {
            read_it_from_str(lineinput, "NKPTS", '=', nkpts);
            if (read_it_from_str(lineinput, "k-points in reciprocal "
            "lattice and")) {
                kvec_wgt.resize(nkpts);
                for (int i=0; i<nkpts; i++)
                    fscanf(ifp, "%lf%lf%lf%lf", &db, &db, &db, &kvec_wgt[i]);
                kvec_wgt /= kvec_wgt.sum();
            }
        }

		if (! read_it_from_str(lineinput, "NBANDsGWLOW"))
			read_it_from_str(lineinput, "NBANDS", '=', nbands);
        read_it_from_str(lineinput, "NIONS", '=', nions);
        if (read_it_from_str(lineinput, "ISPIN", '=', ispin)) {
            int n = nkpts * nbands * ispin;
            band_en.resize(n);
            band_occ.resize(n);
            band_kvec.resize(n);
        }

        int idum;
        if (read_it_from_str(lineinput, "alpha+bet", ':', alpha_plus_beta)) {
            int iev=0;
            for (int isp=0; isp<ispin; isp++) {
                for (int ikp=0; ikp<nkpts; ikp++) {
                    while (fgets(lineinput, 256, ifp) != NULL)
                        if (strstr(lineinput, "band No.") != NULL) break;
                    for (int i=0; i<nbands; i++) {
                        fscanf(ifp, "%d%lf%lf", &idum, &band_en[iev], &band_occ[iev]);
                        band_kvec[iev] = ikp;
                        iev++;
                    }
                }
            }
        }
    }

    fclose(ifp);
}


void PDOS::init_dos()
{
    if (fabs(e_min) + fabs(e_max) < tiny) {
        e_min = band_en.min() - 5 * gaussian_width;
        e_max = band_en.max() + 5 * gaussian_width;
    }
}


void PDOS::read_PDOS_LIST()
{
    if (!PRINT_PDOS && !PRINT_PROFILE)
        return;

    FILE* ifp = NULL;
    if (!open_file("LIST", ifp))
        return;

    char str[256], orbit_name[256], spin_name[256];
    int iion, n = 0;
    while (fscanf(ifp, "%s", str) != EOF) {
        if (str[0] == '#') {                      // comment out strings starting with #
            fgets(str, 256, ifp);
            continue;
        }
        PDOS_LIST item(str);
        pdos_list.push_back(item);
        while (fscanf(ifp, "%d%s%s", &iion, orbit_name, spin_name) == 3) {
            pdos_list[n].ions.push_back(iion);
            if (iion > nions) {
                fprintf(stderr, "ion index %d exceeds the range %d; exit\n", iion, nions);
                exit(-1);
            }
            pdos_list[n].orbitals.push_back(orbit_name);
            pdos_list[n].spins.push_back(spin_name);
        }
        n++;
    }

    fclose(ifp);
}


void PDOS::read_PHASE_LIST()
{
    if (!PRINT_PHASE)
        return;

    FILE* ifp = NULL;
    if (!open_file("LIST", ifp))
        return;

    char str[256], orbit_name[256], spin_name[256];
    int kvec, band, n = 0;
    while (fscanf(ifp, "%s", str) != EOF) {
        if (str[0] == '#') {                      // comment out strings starting with #
            fgets(str, 256, ifp);
            continue;
        }
        PHASE_LIST item(str);
        phase_list.push_back(item);
        while (fscanf(ifp, "%d%d%s%s", &kvec, &band, orbit_name, spin_name) == 4) {
            phase_list[n].kvec.push_back(kvec);
            phase_list[n].band.push_back(band);
            phase_list[n].orbitals.push_back(orbit_name);
            phase_list[n].spins.push_back(spin_name);
        }
        n++;
    }

    fclose(ifp);
}


void PDOS::print_dos(string file_name)
{
    FILE* ofp;
    if ((ofp = fopen(file_name.c_str(), "w")) == NULL) {
        fprintf(stderr, "fail to write %s; exit\n", file_name.c_str());
        exit(-1);
    }

    fprintf(stderr, "DOS is being written on %s\n", file_name.c_str());
    fprintf(ofp, "# gaussian_width = %g\n", gaussian_width);
    for (int i=0; i<dos_energy_mesh.size(); i++)
        fprintf(ofp, "%g %g %g %g\n", dos_energy_mesh[i], full_dos[i],
            occupied_dos[i], accumulated_dos[i]);

    fclose(ofp);
}


void PDOS::print_band(string file_name, string spin_type)
{
    FILE* ofp;
    if ((ofp = fopen(file_name.c_str(), "w")) == NULL) {
        fprintf(stderr, "fail to write %s; exit\n", file_name.c_str());
        exit(-1);
    }

    fprintf(stderr, "BAND is being written on %s\n", file_name.c_str());
    for (int s=0; s<ispin; s++) {
        if (spin_type != "tot" && s != spin_table[spin_type])
            continue;
        int is=s*nkpts*nbands;
        for (int k=0; k<nkpts; k++) {
            int ik=is+k*nbands;
            for (int b=0; b<nbands; b++) {
                int ib=ik+b;
                fprintf(ofp, "%d %g %g %g\n", k%folding_period,
                    band_en[ib], band_occ[ib], kvec_wgt[k]);
                // band_en[ib]+alpha_plus_beta, band_occ[ib], kvec_wgt[k]);
            }
        }
    }

    fclose(ofp);
}


void PDOS::print_BAND()
{
    if (!PRINT_BAND)
        return;

    string file_name = (file_header == "") ? "band" : file_header + ".band";
    print_band(file_name, "tot");
    if (ispin == 2) {
        file_name = (file_header == "") ? "up.band" : file_header + ".up.band";
        print_band(file_name, "up");
        file_name = (file_header == "") ? "dn.band" : file_header + ".dn.band";
        print_band(file_name, "dn");
    }
}


void PDOS::print_band_profile(string file_name, valarray<double>& band_wgt)
{
    FILE* ofp;
    if ((ofp = fopen(file_name.c_str(), "w")) == NULL) {
        fprintf(stderr, "fail to write %s; exit\n", file_name.c_str());
        exit(-1);
    }

    fprintf(stderr, "Band Profile is being written on %s\n", file_name.c_str());
    fprintf(ofp, "# profile cutoff = %g\n", profile_cutoff);
    for (int k=0; k<nkpts; k++) {
        valarray<double> kvec_wgt_for_profile(0., nkpts);
        kvec_wgt_for_profile[k] = 1.;
        calc_dos(band_wgt, band_en, band_occ, band_kvec, kvec_wgt_for_profile);
        valarray<double> peak_dos(full_dos.size());
        peak_dos = full_dos;
        for (int j=0; j<peak_dos.size()-1; j++) {
            if (full_dos[j] < full_dos[j+1]) {
                peak_dos[j+1] += peak_dos[j];
                peak_dos[j] = 0;
            }
        }
        for (int j=peak_dos.size()-1; j>0; j--) {
            if (full_dos[j] < full_dos[j-1]) {
                peak_dos[j-1] += peak_dos[j];
                peak_dos[j] = 0;
            }
        }
        double peak_value;
        for (int j=0; j<dos_energy_mesh.size(); j++) {
            if ((peak_value = peak_dos[j] * interval) > profile_cutoff)
                fprintf(ofp, "%d %g %g\n", k, dos_energy_mesh[j],
                    peak_value);
        }
    }

    fclose(ofp);
}


void PDOS::print_DOS()
{
    if (!PRINT_DOS)
        return;

    init_dos();

    valarray<double> band_wgt(3 - ispin, band_en.size());
    calc_dos(band_wgt, band_en, band_occ, band_kvec, kvec_wgt);
    string file_name = (file_header == "") ? "dos" : file_header + ".dos";
    print_dos(file_name);
    if (ispin == 2) {
        int n = band_en.size() / 2;
        file_name = (file_header == "") ? "up.dos" : file_header + ".up.dos";
        valarray<double> band_occ_tmp(band_occ.size());
        for (int i=0; i<n; i++) {
            band_wgt[i] = 1;
            band_wgt[i+n] = 0;
            band_occ_tmp[i] = band_occ[i];
            band_occ_tmp[i+n] = 0;
        }
        calc_dos(band_wgt, band_en, band_occ_tmp, band_kvec, kvec_wgt);
        print_dos(file_name);
        file_name = (file_header == "") ? "dn.dos" : file_header + ".dn.dos";
        for (int i=0; i<n; i++) {
            band_wgt[i] = 0;
            band_wgt[i+n] = 1;
            band_occ_tmp[i] = 0;
            band_occ_tmp[i+n] = band_occ[i+n];
        }
        calc_dos(band_wgt, band_en, band_occ_tmp, band_kvec, kvec_wgt);
        print_dos(file_name);
    }
}


void PDOS::print_PDOS()
{
    if (!PRINT_PDOS && !PRINT_PROFILE)
        return;

    init_dos();

    if (orbital_wgt.size() == 0 || pdos_list.size() == 0) {
        fprintf(stderr, "PDOS analysis is skipped.\n");
        return;
    }

    int nkb = nkpts * nbands;
    int ndata = (nions == 1) ? 1 : nions + 1;
    int nkbi110 = nkb * ndata * 17;
    for (int i=0; i<pdos_list.size(); i++) {
        fprintf(stderr, "\nPDOS and Band Profile for %s\n",
            pdos_list[i].name.c_str());
        string file_name = (file_header == "") ? "" : file_header + ".";
        file_name += pdos_list[i].name + ".dos";
        valarray<double> band_wgt(0., band_en.size());
        for (int j=0; j<pdos_list[i].ions.size(); j++) {
            int iion = pdos_list[i].ions[j];
            int iorbit = orbital_table[pdos_list[i].orbitals[j]];
            if (ispin == 1)
                pdos_list[i].spins[j] = "tot";
            int isp = spin_table[pdos_list[i].spins[j]];
            int itot = orbital_table["tot"];
            fprintf(stderr, "  atom %d orbit %s spin %s\n", iion,
                pdos_list[i].orbitals[j].c_str(), pdos_list[i].spins[j].c_str());
            for (int is=0; is<ispin; is++) {
                if (isp != 2 && is != isp)
                    continue;
                for (int k=0; k<nkpts; k++) {
                    int ik = is * nkbi110 + k * nbands * ndata * 17;
                    for (int b=0; b<nbands; b++) {
                        int ib = ik + b * ndata * 17;
                        int i1 = ib + (iion - 1) * 17 + iorbit;
                        int i2 = ib + (ndata - 1) * 17 + itot;
                        double wgt = ((orbital_wgt[i2] < 1e-10) ? 0
                            : orbital_wgt[i1] / orbital_wgt[i2]);
                        band_wgt[is*nkb+k*nbands+b] += (3-ispin)*wgt;
                    }
                }
            }
        }
        calc_dos(band_wgt, band_en, band_occ, band_kvec, kvec_wgt);
        if (PRINT_PDOS)
            print_dos(file_name);

        if (PRINT_PROFILE) {
            file_name.replace(file_name.find(".dos"), 5, ".prof");
            print_band_profile(file_name, band_wgt);
        }
    }
}


void PDOS::print_PHASE()
{
    if (!PRINT_PHASE)
        return;

    if (wf_real.size() == 0 || phase_list.size() == 0) {
        fprintf(stderr, "Phase analysis is skipped.\n");
        return;
    }

    int nkb = nkpts * nbands;
    int ndata = nions;
    int nkbi110 = nkb * ndata * 16;

    for (int i=0; i<phase_list.size(); i++) {
        fprintf(stderr, "\nPhase analysis for %s\n", phase_list[i].name.c_str());
        string file_name = (file_header == "") ? "" : file_header + ".";
        file_name += phase_list[i].name + ".wf";
        FILE* ofp;
        if ((ofp = fopen(file_name.c_str(), "w")) == NULL) {
            fprintf(stderr, "fail to write %s; exit\n", file_name.c_str());
            exit(-1);
        }
        fprintf(stderr, "Phase is being written on %s\n", file_name.c_str());
        for (int j=0; j<phase_list[i].kvec.size(); j++) {
            int k = phase_list[i].kvec[j];
            int iband = phase_list[i].band[j];
            int iorbit = orbital_table[phase_list[i].orbitals[j]];
            if (ispin == 1)
                phase_list[i].spins[j] = "tot";
            int isp = spin_table[phase_list[i].spins[j]];
            fprintf(stderr, "kvec %d iband %d orbit %s spin %s\n", k, iband,
                phase_list[i].orbitals[j].c_str(), phase_list[i].spins[j].c_str());
            for (int is=0; is<ispin; is++) {
                if (isp != 2 && is != isp)
                    continue;
                int i1 = (nkpts * is + k - 1) * nbands + iband - 1;
                double ev = band_en[i1];
                double occ = band_occ[i1];
                fprintf(ofp, "# kvec %d iband %d orbit %d isp %d ev %g occ %g\n",
                    k, iband, iorbit, isp, ev, occ);
                fprintf(ofp, "# k = %12.8f%12.8f%12.8f\n",
                    kvec[k-1].pt[0], kvec[k-1].pt[1], kvec[k-1].pt[2]);
                int ik = is * nkbi110 + (k - 1) * nbands * ndata * 16;
                int ib = ik + (iband - 1)* ndata * 16;
                for (int iion=0; iion<nions; iion++) {
                    int l = ib + iion * 16 + iorbit;
                    double rho, theta;
                    complex_conversion(wf_real[l], wf_imag[l], rho, theta);
                    double R, G, B;
                    angle2color(theta, R, G, B);
                    fprintf(ofp, "%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f\n",
                        wf_real[l], wf_imag[l], rho, theta, R, G, B);
                }
            }
        }
        fclose(ofp);
    }
}


void PDOS::complex_conversion(double real, double imag, double& rho, double& theta)
{
    complex<double> c = complex<double>(real, imag);
    rho = abs(c);
    theta = arg(c);
}


void PDOS::angle2color(double angle, double& R, double& G, double& B)
{
    angle -= floor(angle / (2 * M_PI)) * (2 * M_PI);
    double sixty_degree = M_PI / 3;
    if (angle < sixty_degree) {
        R = 1;
        G = (angle - 0 * sixty_degree) / sixty_degree;
        B = 0;
    }
    else if (angle < 2 * sixty_degree) {
        R = (2 * sixty_degree - angle) / sixty_degree;
        G = 1;
        B = 0;
    }
    else if (angle < 3 * sixty_degree) {
        R = 0;
        G = 1;
        B = (angle - 2 * sixty_degree) / sixty_degree;
    }
    else if (angle < 4 * sixty_degree) {
        R = 0;
        G = (4 * sixty_degree - angle) / sixty_degree;
        B = 1;
    }
    else if (angle < 5 * sixty_degree) {
        R = (angle - 4 * sixty_degree) / sixty_degree;
        G = 0;
        B = 1;
    }
    else {
        R = 1;
        G = 0;
        B = (6 * sixty_degree - angle) / sixty_degree;
    }
}


void PDOS::read_PROCAR()
{
    if (!PRINT_PDOS && !PRINT_PROFILE && !PRINT_PHASE)
        return;

    FILE* ifp = NULL;
    if (!open_file("PROCAR", ifp))
        return;

    char *strptr;
    strptr = read_it_from_file(ifp, "PROCAR lm decomposed");
    phase = (strstr(strptr, "phase") != NULL);
    strptr = read_it_from_file(ifp, "# of k-points", ':', nkpts);
    read_it_from_str(strptr, "# of bands", ':', nbands);
    read_it_from_str(strptr, "# of ions", ':', nions);
    int ndata = (nions == 1) ? 1 : nions + 1;     /* +1 for tot */
                                                  /* 17 orbitals, 2 spins */
    orbital_wgt.resize(nkpts * nbands * ndata * 17 * 2);
    kvec.resize(nkpts);
                                                  /* 16 orbitals, 2 spins */
    wf_real.resize(nkpts * nbands * nions * 16 * 2);
                                                  /* 16 orbitals, 2 spins */
    wf_imag.resize(nkpts * nbands * nions * 16 * 2);

    ispin = 1;
    read_orbital_wgt_from_PROCAR(ifp, orbital_wgt, wf_real, wf_imag, nkpts, nbands, nions, ispin);

    if (read_it_from_file(ifp, "# of k-points", ':', nkpts)) {
        fprintf(stderr, "spin-polarized calculation\n");
        ispin = 2;
        read_orbital_wgt_from_PROCAR(ifp, orbital_wgt, wf_real, wf_imag, nkpts, nbands, nions, ispin);
    }

    fclose(ifp);
}


void PDOS::read_orbital_wgt_from_PROCAR(FILE* ifp, valarray<float>& orbital_wgt,
valarray<float>& wf_real, valarray<float>& wf_imag,
int nkpts, int nbands, int nions, int ispin)
{
    int ndata = (nions == 1) ? 1 : nions + 1;
    int indx = nkpts * nbands * ndata * 17 * (ispin - 1);
    int indx1 = nkpts * nbands * nions * 16 * (ispin - 1);
    int indx2 = nkpts * nbands * nions * 16 * (ispin - 1);
    char str[256];
    for (int k=0; k<nkpts; k++) {
        char* strptr;
        strptr = read_it_from_file(ifp, "k-point");
        sscanf(strptr+16, "%lf%lf%lf", &kvec[k].pt[0], &kvec[k].pt[1], &kvec[k].pt[2]);
        for (int j=0; j<nbands; j++) {
            read_it_from_file(ifp, "ion      s     py     pz     px    dxy    dyz ");
            for (int i=0; i<ndata; i++) {
                fscanf(ifp, "%s", &str);
                for (int h=0; h<17; h++)
                    fscanf(ifp, "%f", &orbital_wgt[indx++]);
            }
            if (!phase) continue;
            read_it_from_file(ifp, "ion      s     py     pz     px    dxy    dyz ");
            for (int i=0; i<nions; i++) {
                fscanf(ifp, "%s", &str);
                for (int h=0; h<16; h++)
                    fscanf(ifp, "%f", &wf_real[indx1++]);
                fscanf(ifp, "%s", &str);
                for (int h=0; h<16; h++)
                    fscanf(ifp, "%f", &wf_imag[indx2++]);
            }
        }
    }
}


void PDOS::calc_dos(valarray<double>& band_wgt, valarray<double>& band_en,
valarray<double>& band_occ, valarray<int>& band_kvec,
valarray<double>& kvec_wgt)
{
    int imin = (int)floor((e_min + tiny) / interval);
    int imax = (int)ceil((e_max - tiny)/ interval);
    e_min = imin * interval;
    e_max = imax * interval;

    int n = imax - imin + 1;

    dos_energy_mesh.resize(n);
    for (int i=0; i<n; i++)
        dos_energy_mesh[i] = (imin + i) * interval;
    accumulated_dos.resize(n);
    accumulated_dos = 0;
    valarray<double> accumulated_occupied_dos(n); // temporary
    accumulated_occupied_dos = 0;

    for (int i=0; i<band_en.size(); i++) {
        int immin = (int)floor((band_en[i] - 5 * gaussian_width) / interval);
        int immax = (int)ceil((band_en[i] + 5 * gaussian_width) / interval);
        immin = immin < imin ? imin : (immin > imax ? imax : immin);
        immax = immax < imin ? imin : (immax > imax ? imax : immax);
        double wgt = band_wgt[i] * kvec_wgt[band_kvec[i]];
        double occupied_wgt = band_occ[i] * kvec_wgt[band_kvec[i]];
        for (int j=immin; j<=immax; j++) {
            double x = 1 - 0.5 * erfcc((dos_energy_mesh[j-imin] - band_en[i]
                + 0.5 * interval) /(M_SQRT2 * gaussian_width));
            accumulated_dos[j-imin] += wgt * x;
            accumulated_occupied_dos[j-imin] += occupied_wgt * x;
        }
        for (int j=immax+1; j<=imax; j++) {
            accumulated_dos[j-imin] += wgt;
            accumulated_occupied_dos[j-imin] += occupied_wgt;
        }
    }

    full_dos.resize(n);
    for (int j=1; j<n; j++)
        full_dos[j] = (accumulated_dos[j] - accumulated_dos[j - 1]) / interval;
                                                  // exterapolation
    full_dos[0] = 3*full_dos[1]-3*full_dos[2]+full_dos[3];

    occupied_dos.resize(n);
    for (int j=1; j<n; j++)
        occupied_dos[j] = (accumulated_occupied_dos[j]
            - accumulated_occupied_dos[j - 1]) / interval;
    occupied_dos[0] = 3 * occupied_dos[1] - 3 * occupied_dos[2]
        + occupied_dos[3];                        // exterapolation
}


double PDOS::erfcc(double x)
// Returns the complementary error function erfc(x)
//   erfc(x) = 2 / sqrt(pi) * integral_x^infinity exp(-t^2) dt
// with fractional error everywhere less than 1.2x10^-7
// from Numerical Recipes
{
    double t,z,ans;
    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
        t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
        t*(-0.82215223+t*0.17087277)))))))));
    return x >= 0.0 ? ans : 2.0-ans;
}


template<class T>
bool PDOS::read_it_from_str(char* str, char* to_find, char delimiter, T& value)
{
    char *strptr;
    if ((strptr = strstr(str, to_find)) != NULL) {
        if (delimiter == 0) return true;
        while (*strptr != delimiter && *strptr != '\0') strptr++;
        if (*strptr == '\0') return false;
        istringstream iss(++strptr);
        iss >> value;
        return true;
    }
    return false;
}


template<class T>
char* PDOS::read_it_from_file(FILE* ifp, char* to_find, char delimiter, T& value)
{
    static char lineinput[256];
    while (fgets(lineinput, 256, ifp) != NULL) {
        if (read_it_from_str(lineinput, to_find, delimiter, value))
            return lineinput;
    }
    return NULL;
}
