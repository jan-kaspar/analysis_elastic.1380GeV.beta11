#include <string>
#include <vector>
#include <map>
#include <cmath>

double timestamp0 = 1360537200;

vector<AlignmentSource> alignmentSources;
Analysis anal;
Environment env;

string unsmearing_file;
string unsmearing_object;

string luminosity_data_file;

void Init_base()
{
	// selection of bunches
	keepAllBunches = false;
	bunchMap[9008].push_back(0);
	bunchMap[9009].push_back(0);
	bunchMap[9010].push_back(0);

	// environment settings
	env.InitNominal();
	env.UseMatchedOptics();

	// binning
	// TODO
	anal.t_min = 0.; anal.t_max = 0.5;
	anal.t_min_full = 0.; anal.t_max_full = 0.8;
	
	// approximate (time independent) resolutions - TODO
	anal.si_th_y_1arm = 21E-6 / sqrt(2.);
	anal.si_th_y_1arm_unc = 0.05E-6 / sqrt(2.);

	anal.si_th_y_2arm = anal.si_th_y_1arm / sqrt(2.);
	anal.si_th_y_2arm_unc = 0E-6;

	anal.si_th_x_1arm_L = 10E-6;
	anal.si_th_x_1arm_R = 10E-6;
	anal.si_th_x_1arm_unc = 0E-6;

	anal.si_th_x_2arm = 0E-6;
	anal.si_th_x_2arm_unc = 0E-6;

	// analysis settings
	anal.th_y_lcut = anal.th_y_lcut_L = anal.th_y_lcut_R = 0E-6;
	anal.th_y_hcut = anal.th_y_hcut_L = anal.th_y_hcut_R = 1000E-6;
	
	anal.th_x_lcut = -1.;	
	anal.th_x_hcut = +1.;
	
	anal.use_time_dependent_resolutions = false;

	anal.use_3outof4_efficiency_fits = false;
	anal.use_pileup_efficiency_fits = false;
	anal.inefficiency_3outof4 = 0.;				// diagonal dependent
	anal.inefficiency_shower_near = 0.;
	anal.inefficiency_pile_up = 0.;				// diagonal dependent
	anal.inefficiency_trigger = 0.;

	anal.bckg_corr = 1.;
	
	anal.L_int_eff = 0.;	// mb^-1, diagonal dependent
	
	anal.alignment_t0 = 6000.;		// beginning of the first time-slice
	anal.alignment_ts = 30.*60.;	// time-slice in s
	
	anal.eff_th_y_min = 200E-6; // TODO
	
	anal.t_min_fit = 0.027; // TODO

	// TODO: check
	anal.alignmentYRanges["L_F"] = Analysis::AlignmentYRange(-8., -3.0, 3.6, +8.);
	anal.alignmentYRanges["L_N"] = Analysis::AlignmentYRange(-8., -3.0, 3.4, +8.);
	anal.alignmentYRanges["R_N"] = Analysis::AlignmentYRange(-7., -3.2, 3.0, +7.);
	anal.alignmentYRanges["R_F"] = Analysis::AlignmentYRange(-7., -3.2, 3.0, +7.);

	// TODO
	unsmearing_file = "";	// diagonal dependent
	//unsmearing_object = "cf,<binning>/exp3/corr_final";
	//unsmearing_object = "cf,<binning>/exp3+exp4/corr_final";
	unsmearing_object = "ff";

	// TODO
	luminosity_data_file = "../fill_3549_lumiCalc2.py_V04-02-04_lumibylsXing.csv";
}

//----------------------------------------------------------------------------------------------------

void Init_45b_56t()
{
	// fine alignment settings
	AlignmentSource alSrc;
	alSrc.SetAlignmentA(atConstant);
	alSrc.SetAlignmentB(atConstant);
	alSrc.SetAlignmentC(atConstant);

	// c correction
	//	1st column: from near-far relative alignment
	//	2nd column: from left-right relative alignment
	//	3rd column: from theta*_x vs theta*_y alignment
	alSrc.cnst.a_L_F = 0E-3; alSrc.cnst.b_L_F = 0E-3; alSrc.cnst.c_L_F = +3E-3	+ 110E-3	+ 400E-3;
	alSrc.cnst.a_L_N = 0E-3; alSrc.cnst.b_L_N = 0E-3; alSrc.cnst.c_L_N = -3E-3	+ 110E-3	+ 400E-3;
	alSrc.cnst.a_R_N = 0E-3; alSrc.cnst.b_R_N = 0E-3; alSrc.cnst.c_R_N = -2E-3	+ 110E-3	- 400E-3;
	alSrc.cnst.a_R_F = 0E-3; alSrc.cnst.b_R_F = 0E-3; alSrc.cnst.c_R_F = +2E-3	+ 110E-3	- 400E-3;

	alignmentSources.push_back(alSrc);

	// analysis settings
	anal.cut1_a = 1.; anal.cut1_c = +1E-6; anal.cut1_si = 26E-6;
	anal.cut2_a = 1.; anal.cut2_c = +7E-6; anal.cut2_si = 30E-6;

	anal.cut3_a = +451.; anal.cut3_b = 0.; anal.cut3_si = 0.09;
	anal.cut4_a = -476.; anal.cut4_b = 0.; anal.cut4_si = 0.09;

	anal.cut5_a = -0.01; anal.cut5_b = +0.020; anal.cut5_si = 0.032;
	anal.cut6_a = 0.004; anal.cut6_b = -0.010; anal.cut6_si = 0.032;

	anal.cut7_a = 1095.; anal.cut7_c = +0.0; anal.cut7_si = 0.060;
	anal.cut8_a = 0.; anal.cut8_c = +0.0; anal.cut8_si = 0.125;

	anal.th_y_lcut_L = 145E-6; anal.th_y_lcut_R = 163E-6;
	anal.th_y_lcut = 210E-6;

	// TODO
	//unsmearing_file = "unfolding_fit_45b_56t_old.root";

	anal.inefficiency_3outof4 = 0.0; // TODO
	anal.inefficiency_pile_up = 0.0; // TODO

	anal.L_int_eff = 79.42E3;	// TODO	
}

//----------------------------------------------------------------------------------------------------

void Init_45t_56b()
{
	// fine alignment settings
	AlignmentSource alSrc;
	alSrc.SetAlignmentA(atConstant);
	alSrc.SetAlignmentB(atConstant);
	alSrc.SetAlignmentC(atConstant);

	// c correction
	//	1st column: from near-far relative alignment
	//	2nd column: from left-right relative alignment
	//	3rd column: from theta*_x vs theta*_y alignment
	alSrc.cnst.a_L_F = 0E-3; alSrc.cnst.b_L_F = 0E-3; alSrc.cnst.c_L_F = +100E-3	+ 59E-3	- 0E-3;
	alSrc.cnst.a_L_N = 0E-3; alSrc.cnst.b_L_N = 0E-3; alSrc.cnst.c_L_N = -100E-3	+ 59E-3	- 0E-3;
	alSrc.cnst.a_R_N = 0E-3; alSrc.cnst.b_R_N = 0E-3; alSrc.cnst.c_R_N = -9E-3		+ 59E-3	+ 0E-3;
	alSrc.cnst.a_R_F = 0E-3; alSrc.cnst.b_R_F = 0E-3; alSrc.cnst.c_R_F = +9E-3		+ 59E-3	+ 0E-3;

	alignmentSources.push_back(alSrc);

	// analysis settings
	anal.cut1_a = 1.; anal.cut1_c = 0E-6; anal.cut1_si = 26E-6;
	anal.cut2_a = 1.; anal.cut2_c = -2E-6; anal.cut2_si = 30E-6;

	anal.cut3_a = +426.; anal.cut3_b = 0.; anal.cut3_si = 0.09;
	anal.cut4_a = -461.; anal.cut4_b = 0.; anal.cut4_si = 0.09;

	anal.cut5_a = -0.01; anal.cut5_b = -0.020; anal.cut5_si = 0.035;
	anal.cut6_a = 0.004; anal.cut6_b = +0.020; anal.cut6_si = 0.035;

	anal.cut7_a = 1073.; anal.cut7_c = -0.; anal.cut7_si = 0.060;
	anal.cut8_a = 0.; anal.cut8_c = +0.0; anal.cut8_si = 0.125;

	anal.th_y_lcut_L = 170E-6; anal.th_y_lcut_R = 173E-6;
	anal.th_y_lcut = 210E-6;

	// TODO
	//unsmearing_file = "unfolding_fit_45b_56t_old.root";

	anal.inefficiency_3outof4 = 0.0; // TODO
	anal.inefficiency_pile_up = 0.0; // TODO

	anal.L_int_eff = 79.42E3;	// TODO	
}
