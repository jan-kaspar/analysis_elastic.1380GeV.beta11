#include "TRandom2.h"
#include "TH1D.h"
#include "TRandom2.h"
#include "TVectorD.h"

#include <string>
#include <vector>
#include <cmath>

#include "../Stat.h"
#include "../common_definitions.h"
#include "../common_algorithms.h"

#include "formulae.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

typedef Kinematics (*Func)(const HitData &, const Environment &);

enum RecoQuantity
{
	rqTh_x = 0, rqTh_x_L, rqTh_x_R, rqTh_x_RLdiff, rqTh_x_LFdiff, rqTh_x_RFdiff, rqTh_x_LtoRratio,
	rqVtx_x = 100, rqVtx_x_L, rqVtx_x_R, rqVtx_x_RLdiff,
	rqTh_y = 200, rqTh_y_L, rqTh_y_R, rqTh_y_RLdiff, rqTh_y_LtoRratio,
	rqVtx_y = 300, rqVtx_y_L, rqVtx_y_R, rqVtx_y_RLdiff,
};

enum SimuBits { sbPitch=1, sbBeamDivergence=2, sbMisalignment=4, sbVertex=8, sbOptics=16 };

enum ErrorMode { emAbs, emRel };

//----------------------------------------------------------------------------------------------------

Environment env_nom;

//----------------------------------------------------------------------------------------------------

void AddString(string &s, const string &add)
{
	if (s.empty())
		s = add;
	else
		s += "," + add;
}

//----------------------------------------------------------------------------------------------------

string SimuCodeToString(unsigned int code)
{
	string s;

	if (code & sbPitch) AddString(s, "pitch");
	if (code & sbBeamDivergence) AddString(s, "beamDiv");
	if (code & sbMisalignment) AddString(s, "misalig");
	if (code & sbVertex) AddString(s, "vertex");
	if (code & sbOptics) AddString(s, "optics");

	return s;
}

//----------------------------------------------------------------------------------------------------

Stat st_opt(16);

Stat st_c_th_vtx_L(2), st_c_th_vtx_R(2);
Stat st_d_th_vtx_L(2), st_d_th_vtx_R(2);

Stat st_ta_th_x(1), st_ta_th_y(1);

Stat opt_de_scale_L(2), opt_de_scale_R(2);

void GenerateOpticsPerturbation(Environment &env)
{
	TVectorD de_opt(16);
	env.ApplyRandomOpticsPerturbations(de_opt);
	st_opt.Fill(de_opt);

	//----------

	double D_x_R = + env_nom.L_x_R_N * env_nom.v_x_R_F - env_nom.L_x_R_F * env_nom.v_x_R_N;
	double cth_R = (env.L_x_R_N * env_nom.v_x_R_F - env.L_x_R_F * env_nom.v_x_R_N) / D_x_R;
	double cvtx_R = (env.v_x_R_N * env_nom.v_x_R_F - env.v_x_R_F * env_nom.v_x_R_N) / D_x_R;

	double dth_R = (env_nom.L_x_R_N * env.L_x_R_F - env_nom.L_x_R_F * env.L_x_R_N) / D_x_R;
	double dvtx_R = (env_nom.L_x_R_N * env.v_x_R_F - env_nom.L_x_R_F * env.v_x_R_N) / D_x_R;

	st_c_th_vtx_R.Fill(cth_R, cvtx_R);
	st_d_th_vtx_R.Fill(dth_R, dvtx_R);

	double D_x_L = + env_nom.L_x_L_N * env_nom.v_x_L_F - env_nom.L_x_L_F * env_nom.v_x_L_N;
	double cth_L = (env.L_x_L_N * env_nom.v_x_L_F - env.L_x_L_F * env_nom.v_x_L_N) / D_x_L;
	double cvtx_L = (env.v_x_L_N * env_nom.v_x_L_F - env.v_x_L_F * env_nom.v_x_L_N) / D_x_L;

	double dth_L = (env_nom.L_x_L_N * env.L_x_L_F - env_nom.L_x_L_F * env.L_x_L_N) / D_x_L;
	double dvtx_L = (env_nom.L_x_L_N * env.v_x_L_F - env_nom.L_x_L_F * env.v_x_L_N) / D_x_L;

	st_c_th_vtx_L.Fill(cth_L, cvtx_L);
	st_d_th_vtx_L.Fill(dth_L, dvtx_L);

	//----------

	double off_ta_x = (env.L_x_R_N - env.L_x_R_F) / env_nom.L_x_R_N;
	double off_th_x = cth_R;

	st_ta_th_x.Fill(off_ta_x - off_th_x);
	
	double off_ta_y = (env.L_y_R_N - env.L_y_R_F) / (env_nom.L_y_R_N - env_nom.L_y_R_F);
	double off_th_y = (env.L_y_R_N / env_nom.L_y_R_N + env.L_y_R_F / env_nom.L_y_R_F) / 2.;

	st_ta_th_y.Fill(off_ta_y - off_th_y);

	//----------

	double de_scale_x_L = cth_L - 1.;
	double de_scale_x_R = cth_R - 1.;

	double de_scale_y_L = (env.L_y_L_N / env_nom.L_y_L_N + env.L_y_L_F / env_nom.L_y_L_F) / 2. - 1.;
	double de_scale_y_R = (env.L_y_R_N / env_nom.L_y_R_N + env.L_y_R_F / env_nom.L_y_R_F) / 2. - 1.;

	opt_de_scale_L.Fill(de_scale_x_L, de_scale_y_L);
	opt_de_scale_R.Fill(de_scale_x_R, de_scale_y_R);
}

//---------------------------------------------------------------------------------------------------

void DoOneTest(RecoQuantity q, Func f_reco, unsigned int b, ErrorMode m)
{
	// central values
	Kinematics k_cv;
	k_cv.th_x = k_cv.th_x_L = k_cv.th_x_R = 20E-6;	// rad
	k_cv.th_y = k_cv.th_y_L = k_cv.th_y_R = 20E-6;	// rad
	k_cv.vtx_x = k_cv.vtx_y = 0.;					// mm, must stay 0 (some formulae neglect vertex -> can create bias)!

	// effect histogram
	TH1D *h_eff = new TH1D("", "", 100, 0., 0.);

	// event loop
	unsigned int N_ev = 1000;
	for (unsigned int ev = 0; ev < N_ev; ev++)
	{
		// true
		Kinematics k_tr = k_cv;
		
		// simulation environment
		Environment env_sim = env_nom;

		// beam divergence + vertex
		Kinematics k_sm = k_tr;
		if (b & sbBeamDivergence)
		{
			k_sm.th_x_L += gRandom->Gaus() * env_sim.si_th_x; k_sm.th_y_L += gRandom->Gaus() * env_sim.si_th_y;
			k_sm.th_x_R += gRandom->Gaus() * env_sim.si_th_x; k_sm.th_y_R += gRandom->Gaus() * env_sim.si_th_y;
		}

		if (b & sbVertex)
		{
			k_sm.vtx_x += gRandom->Gaus() * env_sim.si_vtx_x; k_sm.vtx_y += gRandom->Gaus() * env_sim.si_vtx_y;
		}

		// actual optics
		if (b & sbOptics)
		{
			GenerateOpticsPerturbation(env_sim);
		}

		// proton transport
		HitData h = ProtonTransport(k_sm, env_sim);
		
		// misalignment
		if (b & sbMisalignment) 
		{
			h.x_L_F += gRandom->Gaus() * env_sim.si_de_x; h.y_L_F += gRandom->Gaus() * env_sim.si_de_y;
			h.x_L_N += gRandom->Gaus() * env_sim.si_de_x; h.y_L_N += gRandom->Gaus() * env_sim.si_de_y;
			h.x_R_N += gRandom->Gaus() * env_sim.si_de_x; h.y_R_N += gRandom->Gaus() * env_sim.si_de_y;
			h.x_R_F += gRandom->Gaus() * env_sim.si_de_x; h.y_R_F += gRandom->Gaus() * env_sim.si_de_y;
		}

		// pitch error
		if (b & sbPitch)
		{
			h.x_L_F += gRandom->Gaus() * env_sim.si_de_P_L; h.y_L_F += gRandom->Gaus() * env_sim.si_de_P_L;
			h.x_L_N += gRandom->Gaus() * env_sim.si_de_P_L; h.y_L_N += gRandom->Gaus() * env_sim.si_de_P_L;
			h.x_R_N += gRandom->Gaus() * env_sim.si_de_P_R; h.y_R_N += gRandom->Gaus() * env_sim.si_de_P_R;
			h.x_R_F += gRandom->Gaus() * env_sim.si_de_P_R; h.y_R_F += gRandom->Gaus() * env_sim.si_de_P_R;
		}

		// reconstruction
		Environment env_rec = env_nom;
		Kinematics k_re = f_reco(h, env_rec);

		// extract effect
		double q_tr=0., q_re=0.;
		if (q == rqTh_x_R) { q_tr = k_tr.th_x_R; q_re = k_re.th_x_R; }
		if (q == rqTh_x_L) { q_tr = k_tr.th_x_L; q_re = k_re.th_x_L; }
		if (q == rqTh_x) { q_tr = k_tr.th_x; q_re = k_re.th_x; }
		if (q == rqTh_x_RLdiff) { q_tr = k_tr.th_x_R - k_tr.th_x_L; q_re = k_re.th_x_R - k_re.th_x_L; }
		if (q == rqTh_x_LFdiff) { q_tr = k_tr.th_x_L - k_tr.th_x; q_re = k_re.th_x_L - k_re.th_x; }
		if (q == rqTh_x_RFdiff) { q_tr = k_tr.th_x_R - k_tr.th_x; q_re = k_re.th_x_R - k_re.th_x; }
		if (q == rqTh_x_LtoRratio) { q_tr = k_tr.th_x_L / k_tr.th_x_R; q_re = k_re.th_x_L / k_re.th_x_R; }

		if (q == rqVtx_x_R) { q_tr = k_sm.vtx_x; q_re = k_re.vtx_x_R; }
		if (q == rqVtx_x_L) { q_tr = k_sm.vtx_x; q_re = k_re.vtx_x_L; }
		if (q == rqVtx_x) { q_tr = k_sm.vtx_x; q_re = k_re.vtx_x; }
		if (q == rqVtx_x_RLdiff) { q_tr = 0.; q_re = k_re.vtx_x_R - k_re.vtx_x_L; }

		if (q == rqTh_y_R) { q_tr = k_tr.th_y_R; q_re = k_re.th_y_R; }
		if (q == rqTh_y_L) { q_tr = k_tr.th_y_L; q_re = k_re.th_y_L; }
		if (q == rqTh_y) { q_tr = k_tr.th_y; q_re = k_re.th_y; }
		if (q == rqTh_y_RLdiff) { q_tr = k_tr.th_y_R - k_tr.th_y_L; q_re = k_re.th_y_R - k_re.th_y_L; }
		if (q == rqTh_y_LtoRratio) { q_tr = k_tr.th_y_L / k_tr.th_y_R; q_re = k_re.th_y_L / k_re.th_y_R; }

		if (q == rqVtx_y_R) { q_tr = k_sm.vtx_y; q_re = k_re.vtx_y_R; }
		if (q == rqVtx_y_L) { q_tr = k_sm.vtx_y; q_re = k_re.vtx_y_L; }
		if (q == rqVtx_y) { q_tr = k_sm.vtx_y; q_re = k_re.vtx_y; }
		if (q == rqVtx_y_RLdiff) { q_tr = 0.; q_re = k_re.vtx_y_R - k_re.vtx_y_L; }

		double pq = (m == emAbs) ? q_re - q_tr : q_re / q_tr;
		h_eff->Fill(pq);
	}

	string unit;
	double scale = 1.;
	if (q >= rqTh_x && q < rqVtx_x) { unit = "urad"; scale = 1E6; }
	if (q >= rqVtx_x && q < rqTh_y) { unit = "um"; scale = 1E3; }
	if (q >= rqTh_y) { unit = "urad"; scale = 1E6; }

	if (m == emAbs)
		printf("%24s: mean = %+.3f %s, RMS = %.2f %s\n", SimuCodeToString(b).c_str(),
			h_eff->GetMean() * scale, unit.c_str(), h_eff->GetRMS() * scale, unit.c_str());
	else
		printf("%24s: mean = %+.3f, RMS = %.2f %%\n", SimuCodeToString(b).c_str(),
			h_eff->GetMean(), h_eff->GetRMS() * 100.);

	delete h_eff;
}

//----------------------------------------------------------------------------------------------------

void Test(const string &label, RecoQuantity q, Func f)
{
	printf("\n>> %s\n", label.c_str());

	// TODO: uncomment
	/*
	DoOneTest(q, f, sbPitch, emAbs);
	DoOneTest(q, f, sbBeamDivergence, emAbs);
	DoOneTest(q, f, sbVertex, emAbs);
	DoOneTest(q, f, sbPitch | sbBeamDivergence | sbVertex, emAbs);
	DoOneTest(q, f, sbMisalignment, emAbs);
	*/
	DoOneTest(q, f, sbOptics, emRel);
}

//----------------------------------------------------------------------------------------------------

int main()
{
	gRandom->SetSeed(1);

	// nominal environment
	env_nom.InitNominal();
	
	//env_nom.ApplyOpticsCorrections();
	env_nom.UseMatchedOptics();

	// TODO
	//env_nom.si_de_P_L = 12.2E-3; env_nom.si_de_P_R = 12.2E-3;	// mm	(45 bottom - 56 top)
	//env_nom.si_de_P_L = 11.8E-3; env_nom.si_de_P_R = 10.6E-3;	// mm	(45 top - 56 bottom)
	env_nom.si_de_P_L = 12E-3; env_nom.si_de_P_R = 12E-3;	// mm
	env_nom.Print();

	// run tests
	printf("\n");
	printf("==================== theta_x, single pot ====================\n");
	Test("theta_x,1_arm,one_pot_hit_LF:th_x_L", rqTh_x, theta_x_1_arm_one_pot_hit_LF);
	Test("theta_x,1_arm,one_pot_hit_LN:th_x_L", rqTh_x, theta_x_1_arm_one_pot_hit_LN);
	Test("theta_x,1_arm,one_pot_hit_RN:th_x_L", rqTh_x, theta_x_1_arm_one_pot_hit_RN);
	Test("theta_x,1_arm,one_pot_hit_RF:th_x_L", rqTh_x, theta_x_1_arm_one_pot_hit_RF);

	printf("\n");
	printf("==================== theta_x, single arm ====================\n");
	Test("theta_x,1_arm,hit:th_x_L", rqTh_x_L, theta_x_1_arm_hit);
	Test("theta_x,1_arm,angle:th_x_L", rqTh_x_L, theta_x_1_arm_angle);
	Test("theta_x,1_arm,regr:th_x_L", rqTh_x_L, theta_x_1_arm_regr);
	
	printf("\n");
	printf("==================== theta_x, double arm ====================\n");
	Test("theta_x,2_arm,LRavg_hit:th_x", rqTh_x, theta_x_2_arm_LRavg_hit);
	Test("theta_x,2_arm,LRavg_angle:th_x", rqTh_x, theta_x_2_arm_LRavg_angle);
	Test("theta_x,2_arm,LRavg_regr:th_x", rqTh_x, theta_x_2_arm_LRavg_regr);
	Test("theta_x,2_arm,regr:th_x", rqTh_x, theta_x_2_arm_regr);


	printf("\n");
	printf("==================== vtx_x, single arm ====================\n");
	Test("vtx_x,1_arm,regr:vtx_x_L", rqVtx_x_L, vtx_x_1_arm_regr);


	printf("\n");
	printf("==================== vtx_x, double arm ====================\n");
	Test("vtx_x,2_arm,LRavg_regr:vtx_x", rqVtx_x, vtx_x_2_arm_LRavg_regr);
	Test("vtx_x,2_arm,regr:vtx_x", rqVtx_x, vtx_x_2_arm_regr);

	
	printf("\n");
	printf("==================== theta_y, single arm ====================\n");
	Test("theta_y,1_arm,hit:th_y_L", rqTh_y_L, theta_y_1_arm_hit);
	Test("theta_y,1_arm,angle:th_y_L", rqTh_y_L, theta_y_1_arm_angle);
	Test("theta_y,1_arm,regr:th_y_L", rqTh_y_L, theta_y_1_arm_regr);
	
	printf("\n");
	printf("==================== theta_y, double arm ====================\n");
	Test("theta_y,2_arm,LRavg_hit:th_y", rqTh_y, theta_y_2_arm_LRavg_hit);
	Test("theta_y,2_arm,LRavg_angle:th_y", rqTh_y, theta_y_2_arm_LRavg_angle);
	Test("theta_y,2_arm,LRavg_regr:th_y", rqTh_y, theta_y_2_arm_LRavg_regr);
	Test("theta_y,2_arm,regr:th_y", rqTh_y, theta_y_2_arm_regr);

#if 0
	printf("\n");
	printf("============== theta_y, right-left difference ===============\n");
	Test("theta_y_one_arm_hit, th_y_RLdiff", rqTh_y_RLdiff, theta_y_one_arm_hit);
	Test("theta_y_one_arm_angle, th_y_RLdiff", rqTh_y_RLdiff, theta_y_one_arm_angle);
	Test("theta_y_one_arm_regr, th_y_RLdiff", rqTh_y_RLdiff, theta_y_one_arm_regr);
	
	printf("\n");
	printf("================= theta_y, left/right ratio =================\n");
	Test("theta_y_one_arm_hit, th_y_LtoRratio", rqTh_y_LtoRratio, theta_y_one_arm_hit);
	//Test("theta_y_one_arm_angle, th_y_LtoRratio", rqTh_y_LtoRratio, theta_y_one_arm_angle);
	//Test("theta_y_one_arm_regr, th_y_LtoRratio", rqTh_y_LtoRratio, theta_y_one_arm_regr);
#endif	


	printf("\n");
	printf("==================== vtx_y, single arm ====================\n");
	Test("vtx_y,1_arm,regr:vtx_y_L", rqVtx_y_L, vtx_y_1_arm_regr);

	printf("\n");
	printf("==================== vtx_y, double arm ====================\n");
	Test("vtx_y,2_arm,LRavg_regr:vtx_y", rqVtx_y, vtx_y_2_arm_LRavg_regr);
	Test("vtx_y,2_arm,regr:vtx_y", rqVtx_y, vtx_y_2_arm_regr);
	

	// print statistics
	printf("\n");
	printf("==================== statistics ====================\n");
	printf("\n* optics perturbations\n");
	st_opt.PrintMeanAndRMS();
	st_opt.PrintCorrelation();

	printf("\n");

	printf("* st_c_th_vtx_L\n");
	st_c_th_vtx_L.PrintMeanAndRMS();
	st_c_th_vtx_L.PrintCorrelation();

	printf("* st_c_th_vtx_R\n");
	st_c_th_vtx_R.PrintMeanAndRMS();
	st_c_th_vtx_R.PrintCorrelation();
	
	printf("\n");

	printf("* st_d_th_vtx_L\n");
	st_d_th_vtx_L.PrintMeanAndRMS();
	st_d_th_vtx_L.PrintCorrelation();

	printf("* st_d_th_vtx_R\n");
	st_d_th_vtx_R.PrintMeanAndRMS();
	st_d_th_vtx_R.PrintCorrelation();

	printf("\n");

	printf("* st_ta_th_x\n");
	st_ta_th_x.PrintMeanAndRMS();
	
	printf("* st_ta_th_y\n");
	st_ta_th_y.PrintMeanAndRMS();
	
	printf("\n");

	printf("* opt_de_scale_L\n");
	opt_de_scale_L.PrintMeanAndRMS();
	opt_de_scale_L.PrintCorrelation();

	printf("* opt_de_scale_R\n");
	opt_de_scale_R.PrintMeanAndRMS();
	opt_de_scale_R.PrintCorrelation();
}
