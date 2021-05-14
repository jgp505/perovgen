#ifndef _PDOS_H_
#define _PDOS_H_

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


const double tiny = 1e-6;  // a tiny number


class PDOS_LIST {
	public:
		PDOS_LIST(string str) { name = str; }
		string name;
		vector<int> ions;
		vector<string> orbitals;
		vector<string> spins;
};


class PHASE_LIST {
	public:
		PHASE_LIST(string str) { name = str; }
		string name;
		vector<int> kvec;
		vector<int> band;
		vector<string> orbitals;
		vector<string> spins;
};


class POINT {
	public:
		double pt[3];
};


class PDOS {

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
			int i; return read_it_from_str(str, to_find, 0, i); }

		char* read_it_from_file(FILE* ifp, char* to_find) {
			int i; return read_it_from_file(ifp, to_find, 0, i); }
} ;

#endif // _PDOS_H_
