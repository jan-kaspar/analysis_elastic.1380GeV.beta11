#include "TRandom2.h"
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2D.h"

#include <map>
#include <vector>

#include "../../Profile.h"
#include "../../common_definitions.h"
#include "../../common_algorithms.h"
#include "parameters.h"

#include "../classes.h"
#include "../predefined_scenarios.h"
#include "../random_scenarios.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

TGraph* CropDDF(TGraph *fg, double t_min, double t_max)
{
	printf(">> CropDDF\n");
	
	// determine index range
	int i_min = 0, i_max = fg->GetN() - 1;
	double *X = fg->GetX(), *Y = fg->GetY();
	for (int i = 0; i < fg->GetN(); i++)
	{
		if (X[i] < t_min)
			i_min = i;
		if (X[i] > t_max)
		{
			i_max = i;
			break;
		}
	}

	TGraph *g = new TGraph();
	for (int i = i_min; i <= i_max; i++)
	{
		g->SetPoint(g->GetN(), X[i], Y[i]);
	}

	printf("\tpoints selected: %i\n", i_max - i_min + 1);
	printf("\t\tfirst point: idx = %i, |t| = %.4E GeV^2\n", i_min, X[i_min]);
	printf("\t\tlast point: idx = %i, |t| = %.4E GeV^2\n", i_max, X[i_max]);

	return g;
}

//----------------------------------------------------------------------------------------------------

double icdf_sigma_int = 0.;

TGraph* BuildICDF(TGraph *sigma_int, double t_min, double t_max)
{
	printf(">> BuildICDF\n");

	// prepare normalized inverse c.d.f
	int i_min = 0, i_max = sigma_int->GetN() - 1;	// default index range
	for (int i = 0; i < sigma_int->GetN(); i++)
	{
		double x, y;
		sigma_int->GetPoint(i, x, y);
		if (x < t_min)
			i_min = i;
		if (x > t_max)
		{
			i_max = i;
			break;
		}
	}
	
	double p_min = 0., p_max = 0.;
	double t_min_real = 0., t_max_real = 0.;
	sigma_int->GetPoint(i_min, t_min_real, p_min);
	sigma_int->GetPoint(i_max, t_max_real, p_max);

	printf("\tpoints loaded: %i\n", i_max - i_min + 1);
	printf("\t\tfirst point: idx = %i, |t| = %.2E GeV^2, p = %.2E mb\n", i_min, t_min_real, p_min);
	printf("\t\tlast point: idx = %i, |t| = %.2E GeV^2, p = %.2E mb\n", i_max, t_max_real, p_max);
	printf("\t|t| range: %.2E to %.2E GeV^2\n", t_min_real, t_max_real);
	icdf_sigma_int = p_max - p_min;
	printf("\tcorresponding cross-section: %.2E mb\n", icdf_sigma_int);

	if (i_min >= i_max)
		throw "BuildICDF: i_min >= i_max.";

	TGraph *icdf = new TGraph();
	for (int i = i_min; i <= i_max; i++) {
		double x, y;
		sigma_int->GetPoint(i, x, y);
		icdf->SetPoint(icdf->GetN(), (y - p_min) / (p_max - p_min), x);
	}

	return icdf;
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: mc <option> <option> ...\n");
	printf("OPTIONS:\n");
	printf("\t-generator-mode XXX\tgenerator mode: icdf or weights\n");
	printf("\t-model-file XXX\t\tROOT file with model PDF\n");
	printf("\t-model-object XXX\tmodel PDF object name\n");
	printf("\t-scenario-type XXX\tscenario type: random or predefined\n");
	printf("\t-scenario-label XXX\tscenarion name\n");
	printf("\t-seed XXX\t\trandom seed\n");
	printf("\t-events XXX\t\tnumber of events to simulate\n");
	printf("\t-binning XXX\t\tbinning: ub (uniform), eb (exponential), ob (optimised)\n");
	printf("\t-output XXX\t\toutput file name\n");
	printf("\t-diff-plots\t\tbuild difference plots\n");
	printf("\t-debug\t\t\tturn on debug output\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	// defaults
	enum { gmNone, gmICDF, gmWeights } generator_mode = gmWeights;
	string model_file = "input_distributions/t-distributions-4000GeV.root";
	string model_object = "full range/ppp3/PH/differential cross-section";
	enum { stNone, stPredefined, stRandom } scenario_type = stNone;
	string scenario_label = "";
	unsigned int seed = 1;
	unsigned long N_ev = 1E8;
	string binning_type = "ob";
	string outputFileName = "";
	bool diffPlots = false;
	bool debug = false;

	// process command line
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-generator-mode") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				generator_mode = gmNone;
				if (strcmp(argv[i], "icdf") == 0) generator_mode = gmICDF;
				if (strcmp(argv[i], "weights") == 0) generator_mode = gmWeights;
				if (generator_mode == gmNone) {
					printf("ERROR: unknown scenario type `%s'.\n", argv[i]);
					PrintUsage();
					return 2;
				}
			}
			continue;
		}

		if (strcmp(argv[i], "-model-file") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				model_file = argv[i];
			}
			continue;
		}
		
		if (strcmp(argv[i], "-model-object") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				model_object = argv[i];
			}
			continue;
		}

		if (strcmp(argv[i], "-scenario-type") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				if (strcmp(argv[i], "predefined") == 0) scenario_type = stPredefined;
				if (strcmp(argv[i], "random") == 0) scenario_type = stRandom;
				if (scenario_type == stNone)
				{
					printf("ERROR: unknown scenario type `%s'.\n", argv[i]);
					PrintUsage();
					return 2;
				}
			}
			continue;
		}
		
		if (strcmp(argv[i], "-scenario-label") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				scenario_label = argv[i];
			}
			continue;
		}
		
		if (strcmp(argv[i], "-seed") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				seed = atoi(argv[i]);
			}
			continue;
		}

		if (strcmp(argv[i], "-events") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				N_ev = (unsigned long) atof(argv[i]);
			}
			continue;
		}

		if (strcmp(argv[i], "-binning") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				binning_type = argv[i];
			}
			continue;
		}

		if (strcmp(argv[i], "-output") == 0)
		{
			if (argc-1 > i)
			{
				i++;
				outputFileName = argv[i];
			}
			continue;
		}
		
		if (strcmp(argv[i], "-diff-plots") == 0)
		{
			diffPlots = true;
			continue;
		}
		
		if (strcmp(argv[i], "-debug") == 0)
		{
			debug = true;
			continue;
		}

		printf("ERROR: unknown parameter `%s'.\n", argv[i]);
		PrintUsage();
		return 1;
	}

	printf(">> mc executed with these parameters:\n");
	printf("\tgenerator_mode = %u\n", generator_mode);
	printf("\tmodel_file = %s\n", model_file.c_str());
	printf("\tmodel_object = %s\n", model_object.c_str());
	printf("\tscenario_type = %u\n", scenario_type);
	printf("\tscenario_label = %s\n", scenario_label.c_str());
	printf("\tseed = %u\n", seed);
	printf("\tevents = %lu (%.1E)\n", N_ev, double(N_ev));
	printf("\tbinning = %s\n", binning_type.c_str());
	printf("\toutputFileName = %s\n", outputFileName.c_str());
	printf("\tdiffPlots = %u\n", diffPlots);
	printf("\tdebug = %u\n", debug);
	
	bool inputError = false;
	if (generator_mode == gmNone)
	{
		printf("ERROR: generator_mode not set.\n");
		inputError = true;
	}

	if (scenario_type == stNone)
	{
		printf("ERROR: scenario type not set.\n");
		inputError = true;
	}
	
	if (scenario_label.compare("") == 0)
	{
		printf("ERROR: scenario_label not set.\n");
		inputError = true;
	}
	
	if (outputFileName.compare("") == 0)
	{
		printf("ERROR: outputFileName not set.\n");
		inputError = true;
	}
	
	if (inputError)
	{
		PrintUsage();
		return 5;
	}

	// simulation seed
	gRandom->SetSeed(seed);

	// simulation t-distribution
	TFile *fDist = TFile::Open(model_file.c_str());
	if (!fDist)
	{
		printf("ERROR: model file `%s' can not be opened.\n", model_file.c_str());
		return 10;
	}

	TGraph *g_t_dist = (TGraph *) fDist->Get(model_object.c_str());
	if (!g_t_dist)
	{
		printf("ERROR: model object `%s' can not be loaded.\n", model_object.c_str());
		return 11;
	}

	// load settings from parameters.h
	Init_base();

	double th_y_sign = -1.;
	Init_45t_56b();
	
	// initialize environments
	Environment env_nom = env;		// nominal values from parameters.h

	Environment env_sim = env_nom;	// the true values, used in simulation
	Environment env_rec = env_nom;	// the nominal values

	// reset quantities inadequate for simulation
	anal.cut1_a = 1.000; anal.cut1_c = 0E-6;
	anal.cut2_a = 1.000; anal.cut2_c = 0E-6;

	anal.cut3_a = 0.; anal.cut3_b = 0.;
	anal.cut4_a = 0.; anal.cut4_b = 0.;

	// TODO: true value from optics, once well established
	/*anal.cut5_a = 0.10717;*/ anal.cut5_b = -0E-3;
	/*anal.cut6_a = 0.10727;*/ anal.cut6_b = +0E-3;

	anal.cut7_a = 0.; anal.cut7_c = 0E-3;
	anal.cut8_a = 0.; anal.cut8_c = 0E-3;

	anal.si_th_y_1arm = env_nom.si_th_y;
	anal.si_th_y_2arm = 0.;	// not used anywhere

	double norm_corr = 1.;

	// initialise analysis objects
	anal.BuildCuts();
	anal.n_si = 3.;
	Analysis anal_nom = anal;
	
	Analysis anal_id = anal_nom;	// TODO: describe
	Analysis anal_rec = anal_nom;	// the nominal values except for experimental resolutions, used in reconstruction
	
	// simulation t-range
	double t_min_full = anal.t_min_full;
	double t_max_full = anal.t_max_full;
	
	// TODO: remove
	//t_max_full = 0.4;
	
	double t_range = t_max_full - t_min_full;

	//double t_min = anal.t_min;
	//double t_max = anal.t_max;

	// set the systematic scenario
	if (scenario_type == stPredefined)
	{
		if (SetPredefinedScenario(th_y_sign, env_nom, env_sim, anal_id, anal_rec, scenario_label))
		{
			printf("ERROR: unknown predefined scenario `%s'\n", scenario_label.c_str());
			return 5;
		}
	} else {
		if (GenerateRandomScenario(th_y_sign, env_sim, anal_id, anal_rec, scenario_label))
		{
			printf("ERROR: cannot generate random scenario\n");
			return 6;
		}
	}
	
	// print all environments
	printf("\n\n------------------------------ nominal environment ------------------------------\n");
	env_nom.Print();
	
	printf("\n\n------------------------------ simulation (actual) environment ------------------------------\n");
	env_sim.Print();

	printf("\n\n------------------------------ reconstruction environment ------------------------------\n");
	env_rec.Print();
	
	printf("\n\n------------------------------ ideal analysis settings ------------------------------\n");
	anal_id.Print();
	
	printf("\n\n------------------------------ reconstruction analysis settings ------------------------------\n");
	anal_rec.Print();

	printf("\n\n----------------------------------------------------------------------------------------------\n");
	printf("\n\n");

	// output file
	TFile *outF = new TFile(outputFileName.c_str(), "recreate");
	if (outF->IsZombie())
	{
		printf("ERROR: Can't open output file `%s' for writing.\n", outputFileName.c_str());
		return 1;
	}

	// prepare auxiliary plots
	DiffPlots dp_bsm, dp_idre, dp_re, dp_re_idre;
	CutPlots cp_idre(anal_nom), cp_re(anal_nom);

	// prepare t-histograms
	TH1D *h_t_true, *h_t_idre, *h_t_idre_cuts, *h_t_re, *h_t_re_cuts;

	unsigned int N_bins = 0;
	double *binEdges;
	BuildBinning(anal_nom, binning_type, binEdges, N_bins);

	h_t_true = new TH1D("h_t_true", ";|t|   (GeV^{2})", N_bins, binEdges);
	h_t_idre = new TH1D("h_t_idre", ";|t|   (GeV^{2})", N_bins, binEdges);
	h_t_idre_cuts = new TH1D("h_t_idre_cuts", ";|t|   (GeV^{2})", N_bins, binEdges);
	h_t_re = new TH1D("h_t_re", ";|t|   (GeV^{2})", N_bins, binEdges);
	h_t_re_cuts = new TH1D("h_t_re_cuts", ";|t|   (GeV^{2})", N_bins, binEdges);

	h_t_true->Sumw2();
	h_t_idre->Sumw2();
	h_t_idre_cuts->Sumw2();
	h_t_re->Sumw2();
	h_t_re_cuts->Sumw2();

	// prepare th histograms
	TH2D *h_th_x_th_y_true = new TH2D("h_th_x_th_y_true", ";#theta_{x};#theta_{y}", 300, -600E-6, +600E-6, 300, -600E-6, +600E-6);
	TH2D *h_th_x_th_y_re = new TH2D("h_th_x_th_y_re", ";#theta_{x};#theta_{y}", 300, -600E-6, +600E-6, 300, -600E-6, +600E-6);

	// remove unnecessary points from g_t_dist
	TGraph *ddf = CropDDF(g_t_dist, t_min_full, t_max_full);

	// prepare inverse CDF, if needed
	TGraph *icdf = NULL;
	if (generator_mode == gmICDF)
		icdf = BuildICDF(g_t_dist, t_min_full, t_max_full);

	// reset counters
	unsigned long N_ev_acc = 0, N_ev_acc_cut = 0;

	// simulation loop
	for (unsigned long ev = 0; ev < N_ev; ev++)
	{
		if (debug)
		{
			printf("\n");
			printf("----- EVENT %lu -----\n", ev);
		}

		// ----- generate (true) elastic event -----

		double w = 0.;	// event weight

		Kinematics k_tr;
		k_tr.t = 0.;
		k_tr.phi = th_y_sign * gRandom->Rndm() * M_PI;

		if (generator_mode == gmWeights)
		{
			k_tr.t = gRandom->Rndm() * t_range + t_min_full;
			w = ddf->Eval(k_tr.t);
		}

		if (generator_mode == gmICDF)
		{
			double u = gRandom->Rndm();
			k_tr.t = icdf->Eval(u);

			/*
			double B = 20.;		// GeV^-2
			double t_max = 0.03;
			double be = 1. - exp(-B*t_max);
			k_tr.t = -1./B * log(1.-u*be);
			*/
	
			w = 1.;
		}

		Kinematics k_in = k_tr;
		k_in.TPhiToThetas(env_sim);
		
		Kinematics k_id_in = k_tr;
		k_id_in.TPhiToThetas(env_nom);

		// ----- generate beam smearing and vertex -----

		double vtx_x = gRandom->Gaus() * env_sim.si_vtx_x;
		double vtx_y = gRandom->Gaus() * env_sim.si_vtx_y;
		
		double rg_x_R = gRandom->Gaus();
		double rg_x_L = gRandom->Gaus();
		double rg_y_R = gRandom->Gaus();
		double rg_y_L = gRandom->Gaus();
		
		Kinematics k_sm = k_in;
		k_sm.vtx_x = vtx_x;
		k_sm.vtx_y = vtx_y;
		k_sm.th_x_R += rg_x_R * env_sim.si_th_x;
		k_sm.th_x_L += rg_x_L * env_sim.si_th_x;
		k_sm.th_y_R += rg_y_R * env_sim.si_th_y;
		k_sm.th_y_L += rg_y_L * env_sim.si_th_y;
		k_sm.th_x = (k_sm.th_x_R + k_sm.th_x_L) / 2.;
		k_sm.th_y = (k_sm.th_y_R + k_sm.th_y_L) / 2.;
		k_sm.ThetasToTPhi(env_sim);
	
		Kinematics k_id_sm = k_id_in;
		k_id_sm.vtx_x = vtx_x;
		k_id_sm.vtx_y = vtx_y;
		k_id_sm.th_x_R += rg_x_R * env_nom.si_th_x;
		k_id_sm.th_x_L += rg_x_L * env_nom.si_th_x;
		k_id_sm.th_y_R += rg_y_R * env_nom.si_th_y;
		k_id_sm.th_y_L += rg_y_L * env_nom.si_th_y;
		k_id_sm.th_x = (k_id_sm.th_x_R + k_id_sm.th_x_L) / 2.;
		k_id_sm.th_y = (k_id_sm.th_y_R + k_id_sm.th_y_L) / 2.;
		k_id_sm.ThetasToTPhi(env_nom);

		// ----- discard events outside HW acceptance -----
		/*
		if (th_y_L_bsm < th_y_lcut_hw_L || th_y_L_bsm > th_y_hcut_hw_L || th_y_R_bsm < th_y_lcut_hw_R || th_y_R_bsm > th_y_hcut_hw_R)
			continue;		
		*/

		// ----- calculate hit positions -----
		
		HitData h_id_op = ProtonTransport(k_id_sm, env_nom);
		HitData h_op = ProtonTransport(k_sm, env_sim);

		// ----- apply discretization smearing -----

		// (to make sure that the same error is applied to both nom and act hits)
		HitData pitchError;
		pitchError.x_L_F = gRandom->Gaus() * env_nom.si_de_P_L; pitchError.y_L_F = gRandom->Gaus() * env_nom.si_de_P_L;
    	pitchError.x_L_N = gRandom->Gaus() * env_nom.si_de_P_L; pitchError.y_L_N = gRandom->Gaus() * env_nom.si_de_P_L;
    	pitchError.x_R_N = gRandom->Gaus() * env_nom.si_de_P_R; pitchError.y_R_N = gRandom->Gaus() * env_nom.si_de_P_R;
    	pitchError.x_R_F = gRandom->Gaus() * env_nom.si_de_P_R; pitchError.y_R_F = gRandom->Gaus() * env_nom.si_de_P_R;

		HitData h_pi = h_op; 
   		h_pi += pitchError;
		
		HitData h_id_pi = h_id_op;
   		h_id_pi += pitchError;

		// ----- apply misalignment -----
	
		HitData h_mis;
		h_mis.x_L_F = h_pi.x_L_F + env_sim.tilt_L_F*h_pi.y_L_F + env_sim.de_x_L_F; h_mis.y_L_F = h_pi.y_L_F + env_sim.de_y_L_F;
		h_mis.x_L_N = h_pi.x_L_N + env_sim.tilt_L_N*h_pi.y_L_N + env_sim.de_x_L_N; h_mis.y_L_N = h_pi.y_L_N + env_sim.de_y_L_N;
		h_mis.x_R_N = h_pi.x_R_N + env_sim.tilt_R_N*h_pi.y_R_N + env_sim.de_x_R_N; h_mis.y_R_N = h_pi.y_R_N + env_sim.de_y_R_N;
		h_mis.x_R_F = h_pi.x_R_F + env_sim.tilt_R_F*h_pi.y_R_F + env_sim.de_x_R_F; h_mis.y_R_F = h_pi.y_R_F + env_sim.de_y_R_F;
		
		// ----- reconstruction -----

		Kinematics k_id_re = DoReconstruction(h_id_pi, env_nom);	// idre = ideal reconstruction   = simulation:nom      + reconstruction:nom
		Kinematics k_re = DoReconstruction(h_mis, env_rec);			// re = realistic reconstruction = simulation:act(sim) + reconstruction:rec

		if (debug) {
			printf("th_x => input = %+8.2f, smeared = %+8.2f, re = %+8.2f\n", k_in.th_x*1E6, k_sm.th_x*1E6, k_re.th_x*1E6);
			printf("th_y => input = %+8.2f, smeared = %+8.2f, re = %+8.2f\n", k_in.th_y*1E6, k_sm.th_y*1E6, k_re.th_y*1E6);
			/*
			printf("th_x => input = %+E, smeared = %+E, re = %+E\n", k_in.th_x, k_sm.th_x, k_re.th_x);
			printf("th_y => input = %+E, smeared = %+E, re = %+E\n", k_in.th_y, k_sm.th_y, k_re.th_y);
			printf("t_x  => input = %+E, smeared = %+E, re = %+E\n", k_in.t_x, k_sm.t_x, k_re.t_x);
			printf("t_y  => input = %+E, smeared = %+E, re = %+E\n", k_in.t_y, k_sm.t_y, k_re.t_y);
			printf("t    => input = %+E, smeared = %+E, re = %+E\n", k_in.t, k_sm.t, k_re.t);
			*/
		}
		
		// ----- acceptance correction -----

		double phi_corr_idre=0., div_corr_idre=0.;
		bool skip_idre = CalculateAcceptanceCorrections(th_y_sign, k_id_re, anal_id, phi_corr_idre, div_corr_idre);
		phi_corr_idre /= 2.; // only upper half in simulation
		
		double phi_corr_re=0., div_corr_re=0.;
		bool skip_re = CalculateAcceptanceCorrections(th_y_sign, k_re, anal_rec, phi_corr_re, div_corr_re);
		phi_corr_re /= 2.; // only upper half in simulation

		/*
		double phi_corr_idre=1., div_corr_idre=1.;
		bool skip_idre = false;
		
		double phi_corr_re=1., div_corr_re=1.;
		bool skip_re = false;
		*/
		
		// ----- evaluate cuts -----

		CutData cd_id_re;
		bool cut_pass_idre = anal_id.EvaluateCuts(h_id_pi, k_id_re, cd_id_re);
		if (diffPlots)
			cp_idre.Fill(cd_id_re, anal_id);
		
		CutData cd_re;
		bool cut_pass_re = anal_rec.EvaluateCuts(h_mis, k_re, cd_re);
		if (diffPlots)
			cp_re.Fill(cd_re, anal_rec);

		
		// ----- fill plots -----

		if (diffPlots) {
			dp_bsm.Fill(k_sm, k_in, w);
			dp_idre.Fill(k_id_re, k_id_in, w);
			dp_re.Fill(k_re, k_in, w);
			dp_re_idre.Fill(k_re, k_id_re, w);
		}

		h_t_true->Fill(k_tr.t, w);
		h_th_x_th_y_true->Fill(k_in.th_x, k_in.th_y, w);

		if (!skip_idre)
			h_t_idre->Fill(k_id_re.t, w * phi_corr_idre * div_corr_idre);
		if (!skip_idre && cut_pass_idre)
			h_t_idre_cuts->Fill(k_id_re.t, w * phi_corr_idre * div_corr_idre);
		
		if (!skip_re) {
			h_t_re->Fill(k_re.t, w * phi_corr_re * div_corr_re * norm_corr);
			h_th_x_th_y_re->Fill(k_re.th_x, k_re.th_y, w * phi_corr_re * div_corr_re * norm_corr);
			N_ev_acc++;
		}

		if (!skip_re && cut_pass_re)
		{
			h_t_re_cuts->Fill(k_re.t, w * phi_corr_re * div_corr_re * norm_corr);
			N_ev_acc_cut++;
		}
	}

	printf(">> counters\n");
	printf("\tN_ev = %lu\n", N_ev);
	printf("\tN_ev_acc = %lu\n", N_ev_acc);
	printf("\tN_ev_acc_cut = %lu\n", N_ev_acc_cut);
	
	// normalization
	if (generator_mode == gmICDF)
	{
		h_t_true->Scale(icdf_sigma_int / N_ev, "width");
		h_t_idre->Scale(icdf_sigma_int / N_ev, "width");
		h_t_idre_cuts->Scale(icdf_sigma_int / N_ev, "width");
		h_t_re->Scale(icdf_sigma_int / N_ev, "width");
		h_t_re_cuts->Scale(icdf_sigma_int / N_ev, "width");
	}

	if (generator_mode == gmWeights)
	{
		h_t_true->Scale(t_range / N_ev, "width");
		h_t_idre->Scale(t_range / N_ev, "width");
		h_t_idre_cuts->Scale(t_range / N_ev, "width");
		h_t_re->Scale(t_range / N_ev, "width");
		h_t_re_cuts->Scale(t_range / N_ev, "width");
	}

	// save plots
	if (diffPlots)
	{
		gDirectory = outF->mkdir("bsm: diff plots");
		dp_bsm.Write();
		
		gDirectory = outF->mkdir("idre: diff plots");
		dp_idre.Write();
		
		gDirectory = outF->mkdir("re: diff plots");
		dp_re.Write();
		
		gDirectory = outF->mkdir("re - idre: diff plots");
		dp_re_idre.Write();
		
		gDirectory = outF->mkdir("idre: cut plots");
		cp_idre.Write();
		
		gDirectory = outF->mkdir("re: cut plots");
		cp_re.Write();
	}

	gDirectory = outF;
	h_t_true->SetLineColor(1);    	h_t_true->Write();
	h_t_idre->SetLineColor(2);    	h_t_idre->Write();
	h_t_idre_cuts->SetLineColor(7);	h_t_idre_cuts->Write();
	h_t_re->SetLineColor(4);      	h_t_re->Write();
	h_t_re_cuts->SetLineColor(9); 	h_t_re_cuts->Write();

	TCanvas *c = new TCanvas("h_t cmparison");
	c->SetLogy(1);
	h_t_true->Draw("");
	h_t_idre->Draw("same");
	h_t_idre_cuts->Draw("same");
	h_t_re->Draw("same");
	h_t_re_cuts->Draw("same");
	c->Write();

	h_th_x_th_y_true->Write();
	h_th_x_th_y_re->Write();

	h_th_x_th_y_re->Divide(h_th_x_th_y_true);
	h_th_x_th_y_re->SetName("h_th_x_th_y_ratio");
	h_th_x_th_y_re->Write();

	// derived plots
	TH1D *unsm_corr_true = new TH1D(*h_t_true);
	unsm_corr_true->SetName("unsm_corr_true");
	unsm_corr_true->Divide(h_t_idre);
	unsm_corr_true->SetLineColor(2);
	unsm_corr_true->Write();
	
	TH1D *eff_cuts_idre = new TH1D(*h_t_idre_cuts);
	eff_cuts_idre->SetName("eff_cuts_idre");
	eff_cuts_idre->Divide(h_t_idre);
	eff_cuts_idre->SetLineColor(2);
	eff_cuts_idre->Write();
	
	TH1D *eff_cuts_re = new TH1D(*h_t_re_cuts);
	eff_cuts_re->SetName("eff_cuts_re");
	eff_cuts_re->Divide(h_t_re);
	eff_cuts_re->SetLineColor(2);
	eff_cuts_re->Write();
	
	TH1D *h_t_re_tuc = new TH1D(*h_t_re);	// reconstructed histogram treated with TRUE Unsmearing Correction
	h_t_re_tuc->SetName("h_t_re_tuc");
	h_t_re_tuc->Multiply(unsm_corr_true);
	h_t_re_tuc->SetLineColor(2);
	h_t_re_tuc->Write();
	
	TH1D *eff_systematics = new TH1D(*h_t_re_tuc);
	eff_systematics->SetName("eff_systematics");
	for (int bi = 1; bi <= eff_systematics->GetNbinsX(); bi++)
	{
		double f = h_t_true->GetBinContent(bi);
		f = (f != 0.) ? 1./f : 0.;
		eff_systematics->SetBinContent(bi, h_t_re_tuc->GetBinContent(bi) * f);
		eff_systematics->SetBinError(bi, h_t_re_tuc->GetBinError(bi) * f);
	}
	eff_systematics->SetLineColor(8);
	eff_systematics->Write();

	/*
	// load standard unsmearing correction
	TFile *sucF = new TFile("/afs/cern.ch/work/j/jkaspar/analyses/elastic/4000GeV,beta1000/simulation/standard_unsmearing_correction.root");
	TH1D *unsm_corr_std = (TH1D *) sucF->Get("fit_models_merged/h_mean");

	// reconstructed histogram treated with STANDARD Unsmearing Correction
	gDirectory = outF;
	TH1D *h_t_re_suc = new TH1D(*h_t_re);
	h_t_re_suc->SetName("h_t_re_suc");
	h_t_re_suc->Multiply(unsm_corr_std);
	h_t_re_suc->SetLineColor(6);
	h_t_re_suc->Write();
	*/

	
	delete outF;
	return 0;
}
