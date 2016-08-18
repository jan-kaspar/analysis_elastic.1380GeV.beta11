#include <string>
#include <vector>
#include <map>
#include <cmath>

bool updated_ntuples = true;

double timestamp0 = 1360623600;

AlignmentSource alSrc;
Analysis anal;
Environment env;

void Init_base()
{
	// constant-alignment parameters
	alSrc.cnst.a_L_F = 0E-3; alSrc.cnst.b_L_F = 0E-3; alSrc.cnst.c_L_F = 0E-3;
	alSrc.cnst.a_L_N = 0E-3; alSrc.cnst.b_L_N = 0E-3; alSrc.cnst.c_L_N = 0E-3;
	alSrc.cnst.a_R_N = 0E-3; alSrc.cnst.b_R_N = 0E-3; alSrc.cnst.c_R_N = 0E-3;
	alSrc.cnst.a_R_F = 0E-3; alSrc.cnst.b_R_F = 0E-3; alSrc.cnst.c_R_F = 0E-3;

	// alignment settings
	alSrc.SetAlignmentA(atNone);
	alSrc.SetAlignmentB(atNone);
	alSrc.SetAlignmentC(atNone);
	//alSrc.SetAlignmentB(atTimeDependent, "alignment_fit.root");
	//alSrc.SetAlignmentC(atTimeDependent, "alignment_fit.root");

	// environment settings
	env.InitNominal();

	// binning
	// TODO
	anal.t_min = 0.; anal.t_max = 0.5;
	anal.t_min_full = 0.; anal.t_max_full = 0.8;

	// analysis settings
	anal.th_y_lcut_L = 0E-6; anal.th_y_lcut_R = 0E-6;
	anal.th_y_hcut_L = 1000E-6; anal.th_y_hcut_R = 1000E-6;
	
	anal.th_x_lcut_L = -1000E-6; anal.th_x_lcut_R = -1000E-6;	
	anal.th_x_hcut_L = +1000E-6; anal.th_x_hcut_R = +1000E-6;
	
	// TODO
	anal.si_th_y_1arm = 21E-6;
	anal.si_th_y_2arm = anal.si_th_y_1arm / sqrt(2.);
	
	anal.si_th_x_1arm_L = 10E-6;
	anal.si_th_x_1arm_R = 10E-6;
	anal.si_th_x_2arm = anal.si_th_y_2arm;	// TODO
	
	// TODO
	anal.use_full_norm_corr = true;
	anal.full_norm_corr = 0.;
	anal.ineff_corr = 0.;
	anal.bckg_corr = 1.;
	anal.L_int_eff = 0.;	// mb^-1
	
	// TODO
	anal.alignment_t0 = 6000.;		// beginning of the first time-slice
	anal.alignment_ts = 15.*60.;	// time-slice in s
	
	// TODO
	anal.eff_th_y_min = 15E-6;
}


//----------------------------------------------------------------------------------------------------

void Init_45b_56t()
{
	anal.cut1_c = -32E-6; anal.cut1_si = 26E-6;
	anal.cut2_c = -9.E-6; anal.cut2_si = 30E-6;

	anal.cut34_si = 0.120;	

	anal.cut5_a = -0.01; anal.cut5_b = +0.02; anal.cut5_si = 0.035;
	anal.cut6_a = 0.004; anal.cut6_b = -0.01; anal.cut6_si = 0.035;

	anal.cut7_al = 975.; anal.cut7_c = +0.038; anal.cut7_si = 0.040;

	anal.th_y_lcut_L = 160E-6; anal.th_y_lcut_R = 160E-6;

	anal.full_norm_corr = 1.;	// TODO
	anal.ineff_corr = 1.227;	// TODO
	anal.L_int_eff = 14.45E3;	// TODO	
}

//----------------------------------------------------------------------------------------------------

void Init_45t_56b()
{
	anal.cut1_c = -32E-6; anal.cut1_si = 29E-6;
	anal.cut2_c = +30E-6; anal.cut2_si = 30E-6;

	anal.cut34_si = 0.120;

	anal.cut5_a = -0.01; anal.cut5_b = -0.04; anal.cut5_si = 0.035;
	anal.cut6_a = 0.004; anal.cut6_b = -0.18; anal.cut6_si = 0.035;

	anal.cut7_al = 975.; anal.cut7_c = -0.04; anal.cut7_si = 0.040;

	anal.th_y_lcut_L = 190E-6; anal.th_y_lcut_R = 170E-6;

	anal.full_norm_corr = 1.;	// TODO
	anal.ineff_corr = 1.244;	// TODO
	anal.L_int_eff = 14.45E3;	// TODO
}
