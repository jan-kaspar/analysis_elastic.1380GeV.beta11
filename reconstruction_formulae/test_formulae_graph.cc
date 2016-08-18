#include "TRandom2.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TRandom2.h"
#include "TVectorD.h"
#include "TFile.h"

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

enum SimuBits { sbPitch=1, sbBeamDivergence=2, sbMisalignment=4, sbVertex=8, sbAngle=16, sbOptics=32 };

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
	if (code & sbAngle) AddString(s, "angle");
	if (code & sbOptics) AddString(s, "optics");

	return s;
}

//----------------------------------------------------------------------------------------------------

Stat st_opt(16);

void GenerateOpticsPerturbation(Environment &env)
{
	TVectorD de_opt(16);
	env.ApplyRandomOpticsPerturbations(de_opt);
	st_opt.Fill(de_opt);
}

//---------------------------------------------------------------------------------------------------

void TestOnePoint(RecoQuantity q, Func f_reco, unsigned int b, const Kinematics k_cv,
		double &mean, double &stddev)
{
	// effect histogram
	TH1D *h_eff = new TH1D("", "", 100, 0., 0.);

	// event loop
	unsigned int N_ev = 4000;
	for (unsigned int ev = 0; ev < N_ev; ev++)
	{
		// true
		Kinematics k_tr = k_cv;

		// generate scattering angles
		if (b & sbAngle)
		{
			double B = 20.;
			double si_th_x = sqrt(1. / (2. * B * env_nom.p * env_nom.p)), si_th_y = si_th_x;
			k_tr.th_x = k_tr.th_x_L = k_tr.th_x_R = gRandom->Gaus() * si_th_x;
			k_tr.th_y = k_tr.th_y_L = k_tr.th_y_R = gRandom->Gaus() * si_th_y;
		}

		// simulation environment
		Environment env_sim = env_nom;

		// beam divergence + vertex
		Kinematics k_sm = k_tr;
		if (b & sbBeamDivergence)
		{
			k_sm.th_x_L += gRandom->Gaus() * env_sim.si_th_x;
			k_sm.th_y_L += gRandom->Gaus() * env_sim.si_th_y;
			k_sm.th_x_R += gRandom->Gaus() * env_sim.si_th_x;
			k_sm.th_y_R += gRandom->Gaus() * env_sim.si_th_y;
		}

		if (b & sbVertex)
		{
			k_sm.vtx_x += gRandom->Gaus() * env_sim.si_vtx_x;
			k_sm.vtx_y += gRandom->Gaus() * env_sim.si_vtx_y;
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

		double pq = q_re - q_tr;
		h_eff->Fill(pq);
	}

	string unit;
	double scale = 1.;
	if (q >= rqTh_x) { unit = "urad"; scale = 1E6; }
	if (q >= rqVtx_x) { unit = "um"; scale = 1E3; }
	if (q >= rqTh_y) { unit = "urad"; scale = 1E6; }
	if (q >= rqVtx_y) { unit = "um"; scale = 1E3; }

	mean = h_eff->GetMean() * scale;
	stddev = h_eff->GetRMS() * scale;

	delete h_eff;
}

//---------------------------------------------------------------------------------------------------

void TestOneMode(RecoQuantity q, Func f_reco, unsigned int b)
{
	// make subdirectory
	TDirectory *baseDir = gDirectory;
	gDirectory = baseDir->mkdir(SimuCodeToString(b).c_str());

	TGraph *g_mean = new TGraph(); g_mean->SetName("g_mean");
	TGraph *g_stddev = new TGraph(); g_stddev->SetName("g_stddev");
	
	unsigned int independentQuantity = 0;
	double iq_min=0., iq_max=0., iq_step=0.;
	if (q >= rqTh_x) { independentQuantity = 1; iq_min = -100E-6, iq_max = 400E-6, iq_step = 20E-6; }	// rad
	if (q >= rqVtx_x) { independentQuantity = 2; iq_min = -0.020, iq_max = 0.100, iq_step = 0.005; }	// mm
	if (q >= rqTh_y) { independentQuantity = 3; iq_min = -100E-6, iq_max = 400E-6, iq_step = 20E-6; }	// rad
	if (q >= rqVtx_y) { independentQuantity = 4; iq_min = -0.020, iq_max = 0.100, iq_step = 0.005; }	// mm


	for (double iq = iq_min; iq <= iq_max; iq += iq_step)
	{
		// central values
		Kinematics k_cv;
		k_cv.th_x = k_cv.th_x_L = k_cv.th_x_R = 0.;
		k_cv.th_y = k_cv.th_y_L = k_cv.th_y_R = 0.;
		k_cv.vtx_x = k_cv.vtx_y = 0.;

		if (independentQuantity == 1)
			k_cv.th_x = k_cv.th_x_L = k_cv.th_x_R = iq;
		if (independentQuantity == 2)
			k_cv.vtx_x = iq;
		if (independentQuantity == 3)
			k_cv.th_y = k_cv.th_y_L = k_cv.th_y_R = iq;
		if (independentQuantity == 4)
			k_cv.vtx_y = iq;

		double mean=0., stddev=0.;

		TestOnePoint(q, f_reco, b, k_cv, mean=0., stddev=0.);
		g_mean->SetPoint(g_mean->GetN(), iq, mean);
		g_stddev->SetPoint(g_stddev->GetN(), iq, stddev);
	}
	
	g_mean->Write();
	g_stddev->Write();

	// go back to parent directory
	gDirectory = baseDir;
}

//----------------------------------------------------------------------------------------------------

void Test(const string &label, RecoQuantity q, Func f)
{
	printf("\n>> %s\n", label.c_str());

	TDirectory *baseDir = gDirectory;
	gDirectory = baseDir->mkdir(label.c_str());
	
	unsigned int bComplementary = 0;
	if (q >= rqTh_x) { bComplementary = sbVertex; }
	if (q >= rqVtx_x) { bComplementary = sbAngle; }
	if (q >= rqTh_y) { bComplementary = sbVertex; }
	if (q >= rqVtx_y) { bComplementary = sbAngle; }

	TestOneMode(q, f, sbPitch);
	TestOneMode(q, f, sbBeamDivergence);
	TestOneMode(q, f, bComplementary);
	TestOneMode(q, f, sbPitch | sbBeamDivergence | bComplementary);
	//TestOneMode(q, f, sbMisalignment);
	TestOneMode(q, f, sbOptics);
	TestOneMode(q, f, bComplementary | sbOptics);
	TestOneMode(q, f, sbPitch | sbBeamDivergence | bComplementary | sbOptics);
	
	gDirectory = baseDir;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	gRandom->SetSeed(1);

	TFile *f_out = new TFile("test_formulae_graph.root", "recreate");

	// nominal environment
	env_nom.InitNominal();
	// TODO
	//env_nom.si_de_P_L = 12.2E-3; env_nom.si_de_P_R = 12.2E-3;	// mm	(45 bottom - 56 top)
	//env_nom.si_de_P_L = 11.8E-3; env_nom.si_de_P_R = 10.6E-3;	// mm	(45 top - 56 bottom)
	env_nom.si_de_P_L = 11E-3; env_nom.si_de_P_R = 11E-3;	// mm
	env_nom.Print();

	// run tests
	printf("\n");
/*
	printf("==================== theta_x, single pot ====================\n");
	Test("theta_x,1_arm,one_pot_hit_LF:th_x_L", rqTh_x, theta_x_1_arm_one_pot_hit_LF);
	Test("theta_x,1_arm,one_pot_hit_LN:th_x_L", rqTh_x, theta_x_1_arm_one_pot_hit_LN);
	Test("theta_x,1_arm,one_pot_hit_RN:th_x_L", rqTh_x, theta_x_1_arm_one_pot_hit_RN);
	Test("theta_x,1_arm,one_pot_hit_RF:th_x_L", rqTh_x, theta_x_1_arm_one_pot_hit_RF);
*/

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
	
	
	
	printf("==================== theta_y, single pot ====================\n");
	Test("theta_y,1_arm,one_pot_hit_LF:th_y_L", rqTh_y_L, theta_y_1_arm_one_pot_hit_LF);
	Test("theta_y,1_arm,one_pot_hit_LN:th_y_L", rqTh_y_L, theta_y_1_arm_one_pot_hit_LN);

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

	delete f_out;
}
