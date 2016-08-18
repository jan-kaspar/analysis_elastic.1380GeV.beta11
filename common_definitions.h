#ifndef _common_definitions_h_
#define _common_definitions_h_

#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>

#include "TGraph.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TRandom2.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

enum DiagonalType { dUnknown, d45b_56t, d45t_56b, dCombined, ad45b_56b, ad45t_56t };

DiagonalType diagonal = dUnknown;

double th_y_sign = 0.;

//----------------------------------------------------------------------------------------------------

struct AlignmentData
{
	//	a: xy coupling in rad
	//	b: x shift in mm
	//	c: y shift in mm
	double a_L_F, b_L_F, c_L_F;
	double a_L_N, b_L_N, c_L_N;
	double a_R_N, b_R_N, c_R_N;
	double a_R_F, b_R_F, c_R_F;

	AlignmentData()
	{
		a_L_F = b_L_F = c_L_F = 0.;
		a_L_N = b_L_N = c_L_N = 0.;
		a_R_N = b_R_N = c_R_N = 0.;
		a_R_F = b_R_F = c_R_F = 0.;
	}

	AlignmentData Interpolate(double s_N, double s_F, double s_NH, double s_FH) const
	{
		AlignmentData r;

		r.a_L_F = a_L_N + (a_L_F - a_L_N)/(s_F - s_N) * (s_FH - s_N); r.a_L_N = a_L_N + (a_L_F - a_L_N)/(s_F - s_N) * (s_NH - s_N);
		r.a_R_F = a_R_N + (a_R_F - a_R_N)/(s_F - s_N) * (s_FH - s_N); r.a_R_N = a_R_N + (a_R_F - a_R_N)/(s_F - s_N) * (s_NH - s_N);
                                                                          
		r.b_L_F = b_L_N + (b_L_F - b_L_N)/(s_F - s_N) * (s_FH - s_N); r.b_L_N = b_L_N + (b_L_F - b_L_N)/(s_F - s_N) * (s_NH - s_N);
		r.b_R_F = b_R_N + (b_R_F - b_R_N)/(s_F - s_N) * (s_FH - s_N); r.b_R_N = b_R_N + (b_R_F - b_R_N)/(s_F - s_N) * (s_NH - s_N);
                                                                          
		r.c_L_F = c_L_N + (c_L_F - c_L_N)/(s_F - s_N) * (s_FH - s_N); r.c_L_N = c_L_N + (c_L_F - c_L_N)/(s_F - s_N) * (s_NH - s_N);
		r.c_R_F = c_R_N + (c_R_F - c_R_N)/(s_F - s_N) * (s_FH - s_N); r.c_R_N = c_R_N + (c_R_F - c_R_N)/(s_F - s_N) * (s_NH - s_N);

		return r;
	}
};

//----------------------------------------------------------------------------------------------------

enum AlignmentType { atNone, atConstant, atTimeDependent };

struct AlignmentSource
{
	struct GraphSet {
		TGraph *L_F, *L_N, *R_N, *R_F;
		GraphSet() : L_F(NULL), L_N(NULL), R_N(NULL), R_F(NULL) {}
	} gs_a, gs_b, gs_c;

	AlignmentData cnst;

	AlignmentType type_a, type_b, type_c;
	string src_a, src_b, src_c;

	AlignmentSource() : type_a(atNone), type_b(atNone), type_c(atNone)
	{
	}

	void SetAlignmentA(AlignmentType t, const string &fn = "")
	{
		type_a = t;
		src_a = fn;
	}

	void SetAlignmentB(AlignmentType t, const string &fn = "")
	{
		type_b = t;
		src_b = fn;
	}

	void SetAlignmentC(AlignmentType t, const string &fn = "")
	{
		type_c = t;
		src_c = fn;
	}

	void InitOne(const string label, AlignmentType t, const string &fn, GraphSet &gs, const string &obj)
	{
		printf(">> AlignmentSource::InitOne > alignment `%s': type %u\n", label.c_str(), t);

		if (t == atTimeDependent)
		{
			TFile *alF = new TFile(fn.c_str());

			if (alF->IsZombie())
			{
				printf("\tERROR: cannot open file with alignment graphs.\n");
				delete alF;
				return;
			}
			
			TGraph *g_L_F = (TGraph *) alF->Get(( string("L_F/") + obj).c_str() );
			TGraph *g_L_N = (TGraph *) alF->Get(( string("L_N/") + obj).c_str() );
			TGraph *g_R_N = (TGraph *) alF->Get(( string("R_N/") + obj).c_str() );
			TGraph *g_R_F = (TGraph *) alF->Get(( string("R_F/") + obj).c_str() );

			if (g_L_F && g_L_N && g_R_N && g_R_F)
				printf("\talignment graphs successfully loaded\n");
			else {
				printf("\tERROR: unable to load some alignment graphs\n");
				delete alF;
				return;
			}

			gs.L_F = new TGraph(*g_L_F);
			gs.L_N = new TGraph(*g_L_N);
			gs.R_N = new TGraph(*g_R_N);
			gs.R_F = new TGraph(*g_R_F);

			delete alF;
		}
	}

	void Init(const string &path_prefix = "")
	{
		printf(">> AlignmentSource::Init\n");
		InitOne("a", type_a, path_prefix + src_a, gs_a, "a_fit");
		InitOne("b", type_b, path_prefix + src_b, gs_b, "b_fit");
		InitOne("c", type_c, path_prefix + src_c, gs_c, "c_fit");
	}

	AlignmentData Eval(double timestamp) const
	{
		AlignmentData d;

		if (type_a == atNone) {
			d.a_L_F = 0.; d.a_L_N = 0.; d.a_R_N = 0.; d.a_R_F = 0.;
		}
		if (type_a == atConstant) {
			d.a_L_F = cnst.a_L_F; d.a_L_N = cnst.a_L_N; d.a_R_N = cnst.a_R_N; d.a_R_F = cnst.a_R_F;
		}
		if (type_a == atTimeDependent) {
			d.a_L_F = gs_a.L_F->Eval(timestamp)*1E-3; d.a_L_N = gs_a.L_N->Eval(timestamp)*1E-3; d.a_R_N = gs_a.R_N->Eval(timestamp)*1E-3; d.a_R_F = gs_a.R_F->Eval(timestamp)*1E-3;
		}

		if (type_b == atNone) {
			d.b_L_F = 0.; d.b_L_N = 0.; d.b_R_N = 0.; d.b_R_F = 0.;
		}
		if (type_b == atConstant) {
			d.b_L_F = cnst.b_L_F; d.b_L_N = cnst.b_L_N; d.b_R_N = cnst.b_R_N; d.b_R_F = cnst.b_R_F;
		}
		if (type_b == atTimeDependent) {
			d.b_L_F = gs_b.L_F->Eval(timestamp)*1E-3; d.b_L_N = gs_b.L_N->Eval(timestamp)*1E-3; d.b_R_N = gs_b.R_N->Eval(timestamp)*1E-3; d.b_R_F = gs_b.R_F->Eval(timestamp)*1E-3;
		}
		
		if (type_c == atNone) {
			d.c_L_F = 0.; d.c_L_N = 0.; d.c_R_N = 0.; d.c_R_F = 0.;
		}
		if (type_c == atConstant) {
			d.c_L_F = cnst.c_L_F; d.c_L_N = cnst.c_L_N; d.c_R_N = cnst.c_R_N; d.c_R_F = cnst.c_R_F;
		}
		if (type_c == atTimeDependent) {
			d.c_L_F = gs_c.L_F->Eval(timestamp)*1E-3; d.c_L_N = gs_c.L_N->Eval(timestamp)*1E-3; d.c_R_N = gs_c.R_N->Eval(timestamp)*1E-3; d.c_R_F = gs_c.R_F->Eval(timestamp)*1E-3;
		}

		return d;
	}
};

//----------------------------------------------------------------------------------------------------

struct HitData
{
	// validity flags
	unsigned int v_L_F, v_L_N, v_R_F, v_R_N;

	double x_L_F, y_L_F;	// in mm
	double x_L_N, y_L_N;
	double x_R_N, y_R_N;
	double x_R_F, y_R_F;

	HitData() : x_L_F(0.), y_L_F(0.), x_L_N(0.), y_L_N(0.), x_R_N(0.), y_R_N(0.), x_R_F(0.), y_R_F(0.) {}

	void operator += (const HitData &add)
	{
		x_L_F += add.x_L_F; y_L_F += add.y_L_F;
		x_L_N += add.x_L_N; y_L_N += add.y_L_N;
		x_R_N += add.x_R_N; y_R_N += add.y_R_N;
		x_R_F += add.x_R_F; y_R_F += add.y_R_F;
	}

	HitData ApplyAlignment(const AlignmentData &al) const
	{
		HitData r;

		r.x_L_F = x_L_F - al.a_L_F * y_L_F - al.b_L_F; r.y_L_F = y_L_F - al.c_L_F;
		r.x_L_N = x_L_N - al.a_L_N * y_L_N - al.b_L_N; r.y_L_N = y_L_N - al.c_L_N;
		r.x_R_N = x_R_N - al.a_R_N * y_R_N - al.b_R_N; r.y_R_N = y_R_N - al.c_R_N;
		r.x_R_F = x_R_F - al.a_R_F * y_R_F - al.b_R_F; r.y_R_F = y_R_F - al.c_R_F;

		return r;
	}

	HitData ApplyInterpolatedAlignment(const AlignmentData &a, double sN, double sF) const
	{
		AlignmentData a_int = a.Interpolate(214.628, 220.000, sN, sF);

		return ApplyAlignment(a_int);
	}
};

//----------------------------------------------------------------------------------------------------

struct EventRed
{
	HitData h;	// vertical RPs
	HitData hH;	// horizontal RPs

	unsigned int timestamp;
	unsigned int run_num, bunch_num, event_num, trigger_num;
	unsigned int trigger_bits;
};

//----------------------------------------------------------------------------------------------------

struct Environment
{
	// beam momentum (GeV)
	double p, p_L, p_R;

	// beam momentum uncertainty
	double si_de_p;
	
	// beam divergence
	double si_th_x, si_th_y;		// rad

	// vertex smearing
	double si_vtx_x, si_vtx_y; 		// mm

	// pitch-induced error
	double si_de_P_L, si_de_P_R;	// mm

	// optics
	double dvds_x_L, dvds_x_R;	// mm^-1
	double v_x_L_F, v_x_L_N, v_x_R_N, v_x_R_F;
	
	double dvds_y_L, dvds_y_R;	// mm^-1
	double v_y_L_F, v_y_L_N, v_y_R_N, v_y_R_F;
	
	double dLds_x_L, dLds_x_R;
	double L_x_L_F, L_x_L_N, L_x_R_N, L_x_R_F;	// mm
	
	double dLds_y_L, dLds_y_R;
	double L_y_L_F, L_y_L_N, L_y_R_N, L_y_R_F;	// mm

	// optics perturbation covariance matrices
	// order of elements:
	//		left arm:  v_x_L_N, L_x_L_N, v_y_L_N, L_y_L_N, v_x_L_F, L_x_L_F, v_y_L_F, L_y_L_F
	//		right arm: v_x_R_N, L_x_R_N, v_y_R_N, L_y_R_N, v_x_R_F, L_x_R_F, v_y_R_F, L_y_R_F
	// units: v's in 1, L's in m
	TMatrixDSym opt_cov;

	// optics perturbation generator matrices
	TMatrixD opt_per_gen;

	// alignment uncertainties
	double si_de_x, si_de_y, si_tilt;

	// misalignments (mm)
	double de_x_L_N, de_y_L_N, tilt_L_N;
	double de_x_L_F, de_y_L_F, tilt_L_F;
	double de_x_R_N, de_y_R_N, tilt_R_N;
	double de_x_R_F, de_y_R_F, tilt_R_F;

	Environment() : opt_cov(16), opt_per_gen(16, 16)
	{
	}

	void InitNominal();
	void UseMatchedOptics();

	void PrintOpticsUncertainties() const;

	void Print() const
	{
		printf("p=%E, p_L=%E, p_R=%E\n", p, p_L, p_R);
		printf("\n");
		printf("si_th_x=%E, si_th_y=%E\n", si_th_x, si_th_y);
		printf("si_vtx_x=%E, si_vtx_y=%E\n", si_vtx_x, si_vtx_y);
		printf("si_de_P_L=%E, si_de_P_R=%E\n", si_de_P_L, si_de_P_R);
		printf("\n");
		printf("dvds_x_L=%E, dvds_x_R=%E\n", dvds_x_L, dvds_x_R);
		printf("v_x_L_F=%E, v_x_L_N=%E, v_x_R_N=%E, v_x_R_F=%E\n", v_x_L_F, v_x_L_N, v_x_R_N, v_x_R_F);
		printf("\n");
		printf("dvds_y_L=%E, dvds_y_R=%E\n", dvds_y_L, dvds_y_R);
		printf("v_y_L_F=%E, v_y_L_N=%E, v_y_R_N=%E, v_y_R_F=%E\n", v_y_L_F, v_y_L_N, v_y_R_N, v_y_R_F);
		printf("\n");
		printf("dLds_x_L=%E, dLds_x_R=%E\n", dLds_x_L, dLds_x_R);
		printf("L_x_L_F=%E, L_x_L_N=%E, L_x_R_N=%E, L_x_R_F=%E\n", L_x_L_F, L_x_L_N, L_x_R_N, L_x_R_F);
		printf("\n");
		printf("dLds_y_L=%E, dLds_y_R=%E\n", dLds_y_L, dLds_y_R);
		printf("L_y_L_F=%E, L_y_L_N=%E, L_y_R_N=%E, L_y_R_F=%E\n", L_y_L_F, L_y_L_N, L_y_R_N, L_y_R_F);
		printf("\n");
		printf("si_de_x=%E, si_de_y=%E, si_tilt=%E\n", si_de_x, si_de_y, si_tilt);
		printf("\n");
		printf("de_x_L_N=%E, de_y_L_N=%E, tilt_L_N=%E\n", de_x_L_N, de_y_L_N, tilt_L_N);
		printf("de_x_L_F=%E, de_y_L_F=%E, tilt_L_F=%E\n", de_x_L_F, de_y_L_F, tilt_L_F);
		printf("de_x_R_N=%E, de_y_R_N=%E, tilt_R_N=%E\n", de_x_R_N, de_y_R_N, tilt_R_N);
		printf("de_x_R_F=%E, de_y_R_F=%E, tilt_R_F=%E\n", de_x_R_F, de_y_R_F, tilt_R_F);

		PrintOpticsUncertainties();
	}

	void ApplyRandomOpticsPerturbations(TVectorD &de);

	void ApplyRandomOpticsPerturbations()
	{
		TVectorD de(16);
		ApplyRandomOpticsPerturbations(de);
	}
	
	/// modes counted from 0 to 15
	void ApplyOpticsPerturbationMode(int mode, double coef);
};

//----------------------------------------------------------------------------------------------------

void Environment::ApplyRandomOpticsPerturbations(TVectorD &de)
{
	TVectorD r(16);

	for (unsigned int i = 0; i < 16; i++)
		r(i) = gRandom->Gaus();
	de = opt_per_gen * r;
	
	v_x_L_N += de(0) * 1E0;
	L_x_L_N += de(1) * 1E3;
	v_y_L_N += de(2) * 1E0;
	L_y_L_N += de(3) * 1E3;
	v_x_L_F += de(4) * 1E0;
	L_x_L_F += de(5) * 1E3;
	v_y_L_F += de(6) * 1E0;
	L_y_L_F += de(7) * 1E3;

	v_x_R_N += de(8) * 1E0;
	L_x_R_N += de(9) * 1E3;
	v_y_R_N += de(10) * 1E0;
	L_y_R_N += de(11) * 1E3;
	v_x_R_F += de(12) * 1E0;
	L_x_R_F += de(13) * 1E3;
	v_y_R_F += de(14) * 1E0;
	L_y_R_F += de(15) * 1E3;
}

//----------------------------------------------------------------------------------------------------

void Environment::ApplyOpticsPerturbationMode(int mode, double coef)
{
	printf(">> Environment::ApplyOpticsPerturbationMode\n");

	// prepare correlation matrix
	TMatrixDSym cor(opt_cov);
	TMatrixDSym Sigma(opt_cov);
	for (int i = 0; i < opt_cov.GetNrows(); i++)
		for (int j = 0; j < opt_cov.GetNcols(); j++)
		{
			cor(i, j) /= sqrt( opt_cov(i, i) * opt_cov(j, j) );
			Sigma(i, j) = (i == j) ? sqrt( opt_cov(i, i) ) : 0.;
		}

	// eigen decomposition
	TMatrixDSymEigen eig_decomp(cor);
	TVectorD eig_values(eig_decomp.GetEigenValues());

	// construct mode
	TVectorD vm(opt_cov.GetNrows());
	for (int i = 0; i < opt_cov.GetNrows(); i++)
	{
		double l = eig_values(i);
		double sl = (l > 0.) ? sqrt(l) : 0.;
		vm(i) = (i == mode) ? sl * coef : 0.;
	}

	vm = Sigma * eig_decomp.GetEigenVectors() * vm;

	printf("\tleft arm: mode %u, coefficient %+.3f\n", mode, coef);
	vm.Print();

	v_x_L_N += vm(0) * 1E0;
	L_x_L_N += vm(1) * 1E3;
	v_y_L_N += vm(2) * 1E0;
	L_y_L_N += vm(3) * 1E3;
	v_x_L_F += vm(4) * 1E0;
	L_x_L_F += vm(5) * 1E3;
	v_y_L_F += vm(6) * 1E0;
	L_y_L_F += vm(7) * 1E3;

	v_x_R_N += vm(8) * 1E0;
	L_x_R_N += vm(9) * 1E3;
	v_y_R_N += vm(10) * 1E0;
	L_y_R_N += vm(11) * 1E3;
	v_x_R_F += vm(12) * 1E0;
	L_x_R_F += vm(13) * 1E3;
	v_y_R_F += vm(14) * 1E0;
	L_y_R_F += vm(15) * 1E3;
}

//----------------------------------------------------------------------------------------------------

void Environment::PrintOpticsUncertainties() const
{
	printf("optics uncertainties: left arm\n");
	printf("\tv_x_N: %.4f\n", sqrt(opt_cov(0, 0)));
	printf("\tL_x_N: %.4f m\n", sqrt(opt_cov(1, 1)));
	printf("\tv_y_N: %.4f\n", sqrt(opt_cov(2, 2)));
	printf("\tL_y_N: %.4f m\n", sqrt(opt_cov(3, 3)));
	printf("\tv_x_F: %.4f\n", sqrt(opt_cov(4, 4)));
	printf("\tL_x_F: %.4f m\n", sqrt(opt_cov(5, 5)));
	printf("\tv_y_F: %.4f\n", sqrt(opt_cov(6, 6)));
	printf("\tL_y_F: %.4f m\n", sqrt(opt_cov(7, 7)));

	printf("optics uncertainties: right arm\n");
	printf("\tv_x_N: %.4f\n", sqrt(opt_cov(8, 8)));
	printf("\tL_x_N: %.4f m\n", sqrt(opt_cov(9, 9)));
	printf("\tv_y_N: %.4f\n", sqrt(opt_cov(10, 10)));
	printf("\tL_y_N: %.4f m\n", sqrt(opt_cov(11, 11)));
	printf("\tv_x_F: %.4f\n", sqrt(opt_cov(12, 12)));
	printf("\tL_x_F: %.4f m\n", sqrt(opt_cov(13, 13)));
	printf("\tv_y_F: %.4f\n", sqrt(opt_cov(14, 14)));
	printf("\tL_y_F: %.4f m\n", sqrt(opt_cov(15, 15)));
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Kinematics
{
	double th_x_L_F, th_x_L_N, th_x_R_N, th_x_R_F, th_x_L, th_x_R, th_x;	//	rad
	double th_y_L_F, th_y_L_N, th_y_R_N, th_y_R_F, th_y_L, th_y_R, th_y;	//	rad

	double vtx_x_L_F, vtx_x_L_N, vtx_x_R_N, vtx_x_R_F, vtx_x_L, vtx_x_R, vtx_x;	// in mm
	double vtx_y_L_F, vtx_y_L_N, vtx_y_R_N, vtx_y_R_F, vtx_y_L, vtx_y_R, vtx_y;	// in mm

	double th;				// in rad
	double phi;				// in rad
	double t_x, t_y, t;		// in GeV^2

	void ThetasToTPhi(const Environment &env)
	{
		th = sqrt(th_x*th_x + th_y*th_y);
		t_x = env.p*env.p * th_x * th_x;
		t_y = env.p*env.p * th_y * th_y;
		t = t_x + t_y;
		phi = atan2(th_y, th_x);
	}

	void TPhiToThetas(const Environment &env)
	{
		th = sqrt(t) / env.p;
		th_x_L = th_x_R = th_x = th * cos(phi);
		th_y_L = th_y_R = th_y = th * sin(phi);

		t_x = t * cos(phi) * cos(phi);
		t_y = t * sin(phi) * sin(phi);
	}
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct CutData
{
	double cqa[9];	///< array of quantities qa
	double cqb[9];	///< array of quantities qb
	double cv[9];	///< array of cut quantities v = a*qa + b*qb + c
	bool ct[9];		///< array of flags whether |v| < n_si * si
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Analysis
{
	// binning, |t| in GeV^2
	double t_min, t_max;
	double t_min_full, t_max_full;
	double t_min_fit;

	// elastic selection cuts
	double n_si;

	double cut1_a, cut1_c, cut1_si;
	double cut2_a, cut2_c, cut2_si;
	double cut3_a, cut3_b, cut3_si;
	double cut4_a, cut4_b, cut4_si;
	double cut5_a, cut5_b, cut5_si;
	double cut6_a, cut6_b, cut6_si;
	double cut7_a, cut7_c, cut7_si;
	double cut8_a, cut8_c, cut8_si;

	std::vector< std::pair<double, double> > timeIntervals;

	unsigned int N_cuts;	// number of cuts - indexed from 1!
	string cqaN[9], cqbN[9];
	double cca[9], ccb[9], ccc[9], csi[9];
	std::vector<unsigned int> cuts;	// list of active cuts

	// analysis cuts (rad)
	double th_y_lcut, th_y_lcut_L, th_y_lcut_R;
	double th_y_hcut, th_y_hcut_L, th_y_hcut_R;
	
	double th_x_lcut;
	double th_x_hcut;

	// (un)-smearing parameters
	double si_th_x_1arm_L;
	double si_th_x_1arm_R;
	double si_th_x_1arm_unc;

	double si_th_x_2arm;
	double si_th_x_2arm_unc;

	double si_th_y_1arm;
	double si_th_y_1arm_unc;

	double si_th_y_2arm;
	double si_th_y_2arm_unc;
	
	bool use_time_dependent_resolutions;	// whether to use time-dependent fits

	// efficiency parameters
	bool use_3outof4_efficiency_fits;		// whether to use time-dependent fits of 3-out-of-4 efficiency
	bool use_pileup_efficiency_fits;		// whether to use time-dependent fits of pile-up efficiency

	double inefficiency_3outof4;			// inefficiency from 3-out-of-4 method, used only if use_3outof4_efficiency_fits=false
	double inefficiency_shower_near;		// inefficiency due to shower in near RP
	double inefficiency_pile_up;			// inefficiency due to pile-up, used only if use_pileup_efficiency_fits=false
	double inefficiency_trigger;			// trigger inefficiency

	// normalisation correction to subtract background
	double bckg_corr;

	// effective luminosity (dead-time corrected)
	double L_int_eff;	// mb^-1
	
	// alignment slice parameters
	double alignment_t0;	// beginning of the first time-slice
	double alignment_ts;	// time-slice in s
	
	double eff_th_y_min;

	// y ranges for alignment
	struct AlignmentYRange
	{
		double bot_min, bot_max, top_min, top_max;
		AlignmentYRange(double bmi=0., double bma=0., double tmi=0., double tma=0.) :
			bot_min(bmi), bot_max(bma), top_min(tmi), top_max(tma) {}
	};
	map<std::string, AlignmentYRange> alignmentYRanges;

	void BuildCuts();
	bool EvaluateCuts(const HitData &, const Kinematics &, CutData &) const;

	bool SkipTime(unsigned int timestamp) const
	{
		if (timeIntervals.size() == 0)
			return false;

		bool selected = false;
		for (unsigned int i = 0; i < timeIntervals.size(); i++)
		{
			if (timestamp >= timeIntervals[i].first && timestamp <= timeIntervals[i].second) {
				selected = true;
				break;
			}
		}

		return !selected;
	}

	void Print() const
	{
		printf("t_min=%E, t_max=%E, t_min_full=%E, t_max_full=%E\n", t_min, t_max, t_min_full, t_max_full);
		printf("t_min_fit=%E\n", t_min_fit);

		printf("\n");
		printf("%lu time intervals:\n", timeIntervals.size());
		for (std::vector< std::pair<double, double> >::const_iterator it = timeIntervals.begin(); it != timeIntervals.end(); ++it)
			printf("\tfrom %.1f to %.1f\n", it->first, it->second);

		printf("\n");
		printf("n_si=%E\n", n_si);

		printf("\n");
		printf("cut1_a=%E, cut1_c=%E, cut1_si=%E\n", cut1_a, cut1_c, cut1_si);
		printf("cut2_a=%E, cut2_c=%E, cut2_si=%E\n", cut2_a, cut2_c, cut2_si);
		printf("cut3_a=%E, cut3_b=%E, cut3_si=%E\n", cut3_a, cut3_b, cut3_si);
		printf("cut4_a=%E, cut4_b=%E, cut4_si=%E\n", cut4_a, cut4_b, cut4_si);
		printf("cut5_a=%E, cut5_b=%E, cut5_si=%E\n", cut5_a, cut5_b, cut5_si);
		printf("cut6_a=%E, cut6_b=%E, cut6_si=%E\n", cut6_a, cut6_b, cut6_si);
		printf("cut7_a=%E, cut7_c=%E, cut7_si=%E\n", cut7_a, cut7_c, cut7_si);

		printf("\n");
		printf("cut parameters:\n");
		for (unsigned int i = 1; i <= N_cuts; i++)
		{
			printf("%u| cqaN=%s, cqbN=%s | cca=%E, ccb=%E, ccc=%E, csi=%E\n", i,
				cqaN[i].c_str(), cqbN[i].c_str(), cca[i], ccb[i], ccc[i], csi[i]);
		}

		printf("\n");
		printf("%lu enabled cuts: ", cuts.size());
		for (unsigned int i = 0; i < cuts.size(); i++)
			printf((i == 0) ? "%i" : ", %i", cuts[i]);

		printf("\n");
		printf("th_x_lcut=%E\n", th_x_lcut);
		printf("th_x_hcut=%E\n", th_x_hcut);
		printf("th_y_lcut_L=%E, th_y_lcut_R=%E\n", th_y_lcut_L, th_y_lcut_R);
		printf("th_y_hcut_L=%E, th_y_hcut_R=%E\n", th_y_hcut_L, th_y_hcut_R);

		printf("\n");
		printf("si_th_x_1arm_L=%E, si_th_x_1arm_R=%E, si_th_x_1arm_unc=%E\n", si_th_x_1arm_L, si_th_x_1arm_R, si_th_x_1arm_unc);
		printf("si_th_x_2arm=%E, si_th_x_2arm_unc=%E\n", si_th_x_2arm, si_th_x_2arm_unc);
		printf("si_th_y_1arm=%E, si_th_y_1arm_unc=%E\n", si_th_y_1arm, si_th_y_1arm_unc);
		printf("si_th_y_2arm=%E, si_th_y_2arm_unc=%E\n", si_th_y_2arm, si_th_y_2arm_unc);
		printf("use_time_dependent_resolutions = %i\n", use_time_dependent_resolutions);
	
		printf("\n");
		printf("use_3outof4_efficiency_fits = %i\n", use_3outof4_efficiency_fits);
		printf("use_pileup_efficiency_fits= %i\n", use_pileup_efficiency_fits);
		printf("inefficiency_3outof4 = %.3f\n", inefficiency_3outof4);
		printf("inefficiency_shower_near = %.3f\n", inefficiency_shower_near);
		printf("inefficiency_pile_up = %.3f\n", inefficiency_pile_up);
		printf("inefficiency_trigger = %.3f\n", inefficiency_trigger);
		printf("bckg_corr = %.3f\n", bckg_corr);
		printf("L_int_eff=%E\n", L_int_eff);
	}
};

//----------------------------------------------------------------------------------------------------

#include "common_parameters.h"
#include "common_cuts.h"

#endif
