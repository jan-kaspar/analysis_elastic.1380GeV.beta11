#include <string>
#include <vector>
#include <map>
#include <cmath>

double timestamp0 = 1360796400;

vector<AlignmentSource> alignmentSources;
Analysis anal;
Environment env;

string unsmearing_file;
string unsmearing_object;

string luminosity_data_file;

void Init_base()
{
	// selection of bunches
	keepAllBunches = true;

	// alignment settings
	AlignmentSource alSrc;
	alSrc.SetAlignmentA(atNone);
	alSrc.SetAlignmentB(atNone);
	alSrc.SetAlignmentC(atNone);
	alignmentSources.push_back(alSrc);

	alSrc.cnst.a_L_F = 0E-3; alSrc.cnst.b_L_F = 0E-3; alSrc.cnst.c_L_F = 264E-3;
	alSrc.cnst.a_L_N = 0E-3; alSrc.cnst.b_L_N = 0E-3; alSrc.cnst.c_L_N = 264E-3;
	alSrc.cnst.a_R_N = 0E-3; alSrc.cnst.b_R_N = 0E-3; alSrc.cnst.c_R_N = -41E-3;
	alSrc.cnst.a_R_F = 0E-3; alSrc.cnst.b_R_F = 0E-3; alSrc.cnst.c_R_F = -31E-3;
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
	anal.inefficiency_shower_near = 0.03;
	anal.inefficiency_pile_up = 0.;				// diagonal dependent
	anal.inefficiency_trigger = 0.;

	anal.bckg_corr = 1.;
	
	anal.L_int_eff = 0.;	// mb^-1, diagonal dependent
	
	anal.alignment_t0 = 24500.;		// beginning of the first time-slice
	anal.alignment_ts = 15.*60.;	// time-slice in s
	
	anal.eff_th_y_min = 15E-6; // TODO
	
	anal.t_min_fit = 0.027; // TODO

	// TODO
	anal.alignmentYRanges["L_F"] = Analysis::AlignmentYRange(-26., -8.8, 9.2, 27.);
	anal.alignmentYRanges["L_N"] = Analysis::AlignmentYRange(-23., -7.8, 8.0, 24.);
	anal.alignmentYRanges["R_N"] = Analysis::AlignmentYRange(-23., -7.8, 8.0, 23.);
	anal.alignmentYRanges["R_F"] = Analysis::AlignmentYRange(-26., -8.8, 9.0, 26.);

	// TODO
	unsmearing_file = "";	// diagonal dependent
	//unsmearing_object = "cf,<binning>/exp3/corr_final";
	//unsmearing_object = "cf,<binning>/exp3+exp4/corr_final";
	unsmearing_object = "ff";

	// TODO
	luminosity_data_file = "../fill_3564_lumiCalc2.py_V04-02-08_lumibylsXing.csv";
}

//----------------------------------------------------------------------------------------------------

void Init_45b_56t()
{
	anal.cut1_a = 1.; anal.cut1_c = +3E-6; anal.cut1_si = 26E-6;
	anal.cut2_a = 1.; anal.cut2_c = -14E-6; anal.cut2_si = 30E-6;

	anal.cut3_a = 0.; anal.cut3_b = 0.; anal.cut3_si = 0.17;
	anal.cut4_a = 0.; anal.cut4_b = 0.; anal.cut4_si = 0.17;

	anal.cut5_a = -0.01; anal.cut5_b = +0.04; anal.cut5_si = 0.035;
	anal.cut6_a = 0.004; anal.cut6_b = -0.04; anal.cut6_si = 0.035;

	anal.cut7_a = 881.; anal.cut7_c = +0.0; anal.cut7_si = 0.050;

	anal.th_y_lcut_L = 460E-6; anal.th_y_lcut_R = 485E-6;
	//anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6;
	anal.th_y_lcut = 500E-6;

	// TODO
	//unsmearing_file = "unfolding_fit_45b_56t_old.root";

	anal.inefficiency_3outof4 = 0.0; // TODO
	anal.inefficiency_pile_up = 0.0; // TODO

	anal.L_int_eff = 14.45E3;	// TODO	
}

//----------------------------------------------------------------------------------------------------

void Init_45t_56b()
{
	anal.cut1_a = 1.; anal.cut1_c = +5E-6; anal.cut1_si = 29E-6;
	anal.cut2_a = 1.; anal.cut2_c = +19E-6; anal.cut2_si = 30E-6;

	anal.cut3_a = 0.; anal.cut3_b = 0.; anal.cut3_si = 0.17;
	anal.cut4_a = 0.; anal.cut4_b = 0.; anal.cut4_si = 0.17;

	anal.cut5_a = -0.01; anal.cut5_b = -0.075; anal.cut5_si = 0.035;
	anal.cut6_a = 0.004; anal.cut6_b = -0.13; anal.cut6_si = 0.035;

	anal.cut7_a = 820.; anal.cut7_c = -0.; anal.cut7_si = 0.050;

	anal.th_y_lcut_L = 490E-6; anal.th_y_lcut_R = 495E-6;
	//anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6;
	anal.th_y_lcut = 500E-6;

	// TODO
	//unsmearing_file = "unfolding_fit_45b_56t_old.root";

	anal.inefficiency_3outof4 = 0.0; // TODO
	anal.inefficiency_pile_up = 0.0; // TODO

	anal.L_int_eff = 14.45E3;	// TODO
}
