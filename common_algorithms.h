#ifndef _common_algorithms_h_
#define _common_algorithms_h_

#include "TMath.h"
#include "TF1.h"

#include "common_definitions.h"

//----------------------------------------------------------------------------------------------------

Kinematics DoReconstruction(const HitData &h, const Environment &env)
{
	Kinematics k;

	// single-arm reconstruction
	// theta_x: linear regression
	// theta_y: from hit positions
	// vtx_x: linear regression
	// vtx_y: linear regression
	double D_x_L = - env.L_x_L_N * env.v_x_L_F + env.L_x_L_F * env.v_x_L_N;
	k.th_x_L = (env.v_x_L_F * h.x_L_N - env.v_x_L_N * h.x_L_F) / D_x_L;
	k.vtx_x_L = (env.L_x_L_F * h.x_L_N - env.L_x_L_N * h.x_L_F) / D_x_L;
	
	double D_x_R = env.L_x_R_N * env.v_x_R_F - env.L_x_R_F * env.v_x_R_N;
	k.th_x_R = (env.v_x_R_F * h.x_R_N - env.v_x_R_N * h.x_R_F) / D_x_R;
	k.vtx_x_R = (-env.L_x_R_F * h.x_R_N + env.L_x_R_N * h.x_R_F) / D_x_R;

	k.th_y_L_N = - h.y_L_N / env.L_y_L_N;
  	k.th_y_L_F = - h.y_L_F / env.L_y_L_F;
  	k.th_y_L = (k.th_y_L_N + k.th_y_L_F) / 2.;
  	
	k.th_y_R_N = + h.y_R_N / env.L_y_R_N;
  	k.th_y_R_F = + h.y_R_F / env.L_y_R_F;
  	k.th_y_R = (k.th_y_R_N + k.th_y_R_F) / 2.;
	
	double D_y_L = - env.L_y_L_N * env.v_y_L_F + env.L_y_L_F * env.v_y_L_N;
	//k.th_y_L = (env.v_y_L_F * h.y_L_N - env.v_y_L_N * h.y_L_F) / D_y_L;
	k.vtx_y_L = (env.L_y_L_F * h.y_L_N - env.L_y_L_N * h.y_L_F) / D_y_L;
	
	double D_y_R = env.L_y_R_N * env.v_y_R_F - env.L_y_R_F * env.v_y_R_N;
	//k.th_y_R = (env.v_y_R_F * h.y_R_N - env.v_y_R_N * h.y_R_F) / D_y_R;
	k.vtx_y_R = (-env.L_y_R_F * h.y_R_N + env.L_y_R_N * h.y_R_F) / D_y_R;

	// double-arm reconstruction
	// theta_x: linear regression using all 4 measurements
	// theta_y: L-R average of single-arm hit results
	// vtx_x: L-R average of single-arm regression results
	// vtx_y: L-R average of single-arm regression results
	k.vtx_x = (k.vtx_x_L + k.vtx_x_R) / 2.;
	k.vtx_y = (k.vtx_y_L + k.vtx_y_R) / 2.;
	
	k.th_x = (k.th_x_L + k.th_x_R) / 2.;
	k.th_y = (k.th_y_L + k.th_y_R) / 2.;

	/*
	double SLL_x = + env.L_x_L_F*env.L_x_L_F + env.L_x_L_N*env.L_x_L_N + env.L_x_R_N*env.L_x_R_N + env.L_x_R_F*env.L_x_R_F;
	double SLv_x = - env.L_x_L_F*env.v_x_L_F - env.L_x_L_N*env.v_x_L_N + env.L_x_R_N*env.v_x_R_N + env.L_x_R_F*env.v_x_R_F;
	double Svv_x = + env.v_x_L_F*env.v_x_L_F + env.v_x_L_N*env.v_x_L_N + env.v_x_R_N*env.v_x_R_N + env.v_x_R_F*env.v_x_R_F;
	double D_x = (SLL_x * Svv_x - SLv_x * SLv_x);
	
	double SLh_x = - env.L_x_L_F*h.x_L_F - env.L_x_L_N*h.x_L_N + env.L_x_R_N*h.x_R_N + env.L_x_R_F*h.x_R_F;
	double Svh_x = + env.v_x_L_F*h.x_L_F + env.v_x_L_N*h.x_L_N + env.v_x_R_N*h.x_R_N + env.v_x_R_F*h.x_R_F;
	
	k.th_x = (Svv_x * SLh_x - SLv_x * Svh_x) / D_x;
	*/

	// theta reconstruction
	double th_sq = k.th_x*k.th_x + k.th_y*k.th_y;
	k.th = sqrt(th_sq);
	k.phi = atan2(k.th_y, k.th_x);

	// t reconstruction
	k.t_x = env.p*env.p * k.th_x * k.th_x;
	k.t_y = env.p*env.p * k.th_y * k.th_y;
	k.t = k.t_x + k.t_y;

	return k;
}

//----------------------------------------------------------------------------------------------------

void BuildBinning(const Analysis &anal, const string &type, double* &binEdges, unsigned int &bins,
		const string &prefix_path = "", bool verbose = false)
{
	if (verbose)
		printf(">> BuildBinning(%s)\n", type.c_str());
	
	std::vector<double> be;
	double w;

	// between t_min_full and t_min
	unsigned int N_bins_low = 0;
	w = (anal.t_min - anal.t_min_full) / N_bins_low;
	for (unsigned int i = 0; i < N_bins_low; i++)
		be.push_back(anal.t_min_full + w * i);

	// between t_min and t_max
	unsigned int N_bins_cen = 200;
	
	if (type.compare("ub") == 0)
	{
		w = (anal.t_max - anal.t_min) / N_bins_cen;
		for (unsigned int i = 0; i < N_bins_cen; i++)
			be.push_back(anal.t_min + w * i);
	}
	
	if (type.compare("eb") == 0)
	{
		double B = 3.;
		for (unsigned int bi = 0; bi < N_bins_cen; bi++)
			be.push_back( - log( (1. - double(bi) / N_bins_cen) * exp(-B*anal.t_min) + double(bi) * exp(-B*anal.t_max) / N_bins_cen ) / B );
	}

	if (type.find("ob") == 0)
	{
		// extract parameters
		size_t p1 = type.find("-", 0);
		size_t p2 = type.find("-", p1 + 1);
		size_t p3 = type.find("-", p2 + 1);
		
		double n_smearing_sigmas = atof(type.substr(p1+1, p2-p1-1).c_str());
		string stat_unc_label = type.substr(p2+1, p3-p2-1);
		double bs_max = atof(type.substr(p3+1).c_str());

		// load generators
		TFile *f_in = TFile::Open((prefix_path + "../binning/generators.root").c_str());
		TGraph *g_rms_t = (TGraph *) f_in->Get("g_rms_t");
		TGraph *g_bs_fsu = (TGraph *) f_in->Get( ("g_bs_stat_unc_" + stat_unc_label).c_str() );

		double t = anal.t_min;
		while (t < anal.t_max)
		{
			be.push_back(t);

			double w = max(n_smearing_sigmas * g_rms_t->Eval(t), g_bs_fsu->Eval(t));
			double t_c = t + w/2.;
			w = max(n_smearing_sigmas * g_rms_t->Eval(t_c), g_bs_fsu->Eval(t_c));
			if (w > bs_max)
				w = bs_max;

			t += w;
		}

		delete f_in;
	}

	// between t_max and t_max_full
	unsigned int N_bins_high = 50;
	w = (anal.t_max_full - anal.t_max) / N_bins_high;
	for (unsigned int i = 0; i <= N_bins_high; i++)
		be.push_back(anal.t_max + w * i);

	// return results
	bins = be.size() - 1;
	binEdges = new double[be.size()];
	for (unsigned int i = 0; i < be.size(); i++)
	{
		binEdges[i] = be[i];
		if (verbose)
			printf("\tbi = %4u: %.4E\n", i, binEdges[i]);
	}
}

//----------------------------------------------------------------------------------------------------

bool CalculateAcceptanceCorrections(double th_y_sign,
		const Kinematics &k, const Analysis &anal,
		double &phi_corr, double &div_corr)
{
	// ---------- smearing component ----------

	/*
	if ((k.th_x_L < anal.th_x_lcut_L) || (k.th_x_R < anal.th_x_lcut_R) || (k.th_x_L > anal.th_x_hcut_L) || (k.th_x_R > anal.th_x_hcut_R))
		return true;
	*/

	if ((th_y_sign * k.th_y_L < anal.th_y_lcut_L) || (th_y_sign * k.th_y_R < anal.th_y_lcut_R)
		|| (th_y_sign * k.th_y_L > anal.th_y_hcut_L) || (th_y_sign * k.th_y_R > anal.th_y_hcut_R))
		return true;
	
	/*
	double LB_x_L = anal.th_x_lcut_L - k.th_x, UB_x_L = anal.th_x_hcut_L - k.th_x;
	double LB_x_R = anal.th_x_lcut_R - k.th_x, UB_x_R = anal.th_x_hcut_R - k.th_x;
	double F_x_L = (UB_x_L > LB_x_L) ? ( TMath::Erf(UB_x_L / anal.si_th_x_1arm_L / sqrt(2.)) - TMath::Erf(LB_x_L / anal.si_th_x_1arm_L / sqrt(2.)) ) / 2. : 0.;
	double F_x_R = (UB_x_R > LB_x_R) ? ( TMath::Erf(UB_x_R / anal.si_th_x_1arm_R / sqrt(2.)) - TMath::Erf(LB_x_R / anal.si_th_x_1arm_R / sqrt(2.)) ) / 2. : 0.;
	double F_x = F_x_L * F_x_R;
	*/
	double F_x = 1.;

	double th_y_abs = th_y_sign * k.th_y;

	double UB_y = min(anal.th_y_hcut_R - th_y_abs, th_y_abs - anal.th_y_lcut_L);
	double LB_y = max(anal.th_y_lcut_R - th_y_abs, th_y_abs - anal.th_y_hcut_L);
	double F_y = (UB_y > LB_y) ? ( TMath::Erf(UB_y / anal.si_th_y_1arm) - TMath::Erf(LB_y / anal.si_th_y_1arm) ) / 2. : 0.;

	//printf(">> F_x_L = %E, F_x_R = %E, F_y = %E\n", F_x_L, F_x_R, F_y);

	div_corr = 1./(F_x * F_y);

	// ---------- phi component ----------
	
	// apply safety margins to avoid excessive smearing component
	double th_x_lcut = anal.th_x_lcut;
	double th_x_hcut = anal.th_x_hcut;

	//double th_y_lcut = max(anal.th_y_lcut_L, anal.th_y_lcut_R) + 0.2E-6;
	//double th_y_hcut = min(anal.th_y_hcut_L, anal.th_y_hcut_R) - 1.0E-6;
	double th_y_lcut = anal.th_y_lcut;
	double th_y_hcut = anal.th_y_hcut;

	if (k.th_x <= th_x_lcut || k.th_x >= th_x_hcut || th_y_abs <= th_y_lcut || th_y_abs >= th_y_hcut)
		return true;

	// get all intersections
	set<double> phis;

	if (k.th > th_y_lcut)
	{
		double phi = asin(th_y_lcut / k.th);
		double ta_x = k.th * cos(phi);
		if (th_x_lcut < ta_x && ta_x < th_x_hcut)
			phis.insert(phi);
		if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
			phis.insert(M_PI - phi);
	}
	
	if (k.th > th_y_hcut)
	{
		double phi = asin(th_y_hcut / k.th);
		double ta_x = k.th * cos(phi);
		if (th_x_lcut < ta_x && ta_x < th_x_hcut)
			phis.insert(phi);
		if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
			phis.insert(M_PI - phi);
	}

	if (k.th > fabs(th_x_hcut))
	{
		double phi = acos(fabs(th_x_hcut) / k.th);
		double ta_y = k.th * sin(phi);
		if (th_y_lcut < ta_y && ta_y < th_y_hcut)
			phis.insert(phi);
	}

	if (k.th > fabs(th_x_lcut))
	{
		double phi = acos(fabs(th_x_lcut) / k.th);
		double ta_y = k.th * sin(phi);
		if (th_y_lcut < ta_y && ta_y < th_y_hcut)
			phis.insert(M_PI - phi);
	}

	// the number of intersections must be even
	if ((phis.size() % 2) == 1)
	{
		printf("ERROR: odd number of intersections in acceptance calculation\n");
	}

	// no intersection => no acceptances
	if (phis.size() == 0)
		return true;

	// calculate arc-length in within acceptance
	double phiSum = 0.;
	for (set<double>::iterator it = phis.begin(); it != phis.end(); ++it)
	{
		double phi_start = *it;
		++it;
		double phi_end = *it;

		phiSum += phi_end - phi_start;
	}
	
	phi_corr = 2. * M_PI / phiSum;

	return false;
}

//----------------------------------------------------------------------------------------------------

bool SkipBlindedEvent(unsigned int event_num)
{
	if ((event_num % 2) == 1)
		return true;

	return false;
}

//----------------------------------------------------------------------------------------------------

bool SkipRun(unsigned int /*run*/, unsigned int /* file */, bool /* strict */ = true)
{
	// TODO - why ??
	//if (run == 9008)
	//	return true;

	return false;
}

//----------------------------------------------------------------------------------------------------

// map: run number (8372) --> list of triggered bunches
typedef std::map<unsigned int, std::vector<unsigned int> > BunchMap;

BunchMap bunchMap;
bool keepAllBunches = false;

bool SkipBunch(unsigned int run, unsigned bunch)
{
	if (keepAllBunches)
		return false;

	const std::vector<unsigned int> &bunches = bunchMap[run];
	return (find(bunches.begin(), bunches.end(), bunch) == bunches.end());
}

//----------------------------------------------------------------------------------------------------

// returns the beam for which the bunch is non-colliding
// for colliding bunches returns zero
unsigned int NonCollidingBunch(unsigned int run, unsigned bunch)
{
	if (run >= 9008 && run <= 9010) {
		if (bunch == 894 || bunch == 2773 || bunch == 3000)
			return 1;
		if (bunch == 891 || bunch == 2770 || bunch == 2980)
			return 2;
	}

	return 0;
}

//----------------------------------------------------------------------------------------------------

bool IsZeroBias(unsigned int trigger, unsigned int /* run */, unsigned int /* event */)
{
	return ((trigger & 512) != 0);
}

//----------------------------------------------------------------------------------------------------

HitData ProtonTransport(const Kinematics &k, const Environment &env)
{
	HitData h;
	h.x_L_F = -env.L_x_L_F*k.th_x_L + env.v_x_L_F*k.vtx_x; h.y_L_F = -env.L_y_L_F*k.th_y_L + env.v_y_L_F*k.vtx_y;
	h.x_L_N = -env.L_x_L_N*k.th_x_L + env.v_x_L_N*k.vtx_x; h.y_L_N = -env.L_y_L_N*k.th_y_L + env.v_y_L_N*k.vtx_y;
	h.x_R_N = +env.L_x_R_N*k.th_x_R + env.v_x_R_N*k.vtx_x; h.y_R_N = +env.L_y_R_N*k.th_y_R + env.v_y_R_N*k.vtx_y;
	h.x_R_F = +env.L_x_R_F*k.th_x_R + env.v_x_R_F*k.vtx_x; h.y_R_F = +env.L_y_R_F*k.th_y_R + env.v_y_R_F*k.vtx_y;

	return h;
}

#endif
