#include <cstdio>
#include <cstring>
#include <string>

#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TH1D.h"
#include "TRandom3.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

double df_t(double t)
{
	t *= 1.06;	// in order to get si_th_x_L and si_th_x_R (after additional smearing) 122.9 urad
	return exp(-1.46462e+01 * t -5.92124e+00 *t*t);
}

//----------------------------------------------------------------------------------------------------

TSpline* BuildICDF(double t_min, double t_max)
{
	// make graph of CDF
	TGraph *g = new TGraph();
	double cdf_min = -1.;
	double cdf_max = -1.;
	double t_prev = -1.;
	double df_prev = -1.;
	double cdf = 0.;
	for (double t = t_min; ; t += 1E-3)
	{
		double df = df_t(t);

		if (t_prev > 0.)
			cdf += (df + df_prev) / 2. / (t - t_prev);

		g->SetPoint(g->GetN(), cdf, t);

		if (cdf < cdf_min || cdf_min < 0.)
			cdf_min = cdf;

		if (cdf > cdf_max || cdf_max < 0.)
			cdf_max = cdf;

		if (t >= t_max)
			break;

		t_prev = t;
		df_prev = df;
	}

	// normalise CDF
	for (int i = 0; i < g->GetN(); i++)
	{
		double cdf, t;
		g->GetPoint(i, cdf, t);
		g->SetPoint(i, (cdf - cdf_min) / (cdf_max - cdf_min), t);

		//printf("%i, %.30E, %E\n", i, (cdf - cdf_min) / (cdf_max - cdf_min), v);
	}

	// transform graph to spline
	TSpline3 *s = new TSpline3("", g->GetX(), g->GetY(), g->GetN());
	
	delete g;
	return s;
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: simulation <option> <option> ...\n");
	printf("OPTIONS:\n");
	printf("    --diagonal      \n");
	printf("    --outputDir     \n");
	printf("    --de_th_y       \n");
	printf("    --nEv           \n");
	printf("    --seed          \n");
	printf("    --nonExponential\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string diagonal = "";
	string outputDir = "./";
	unsigned int N_ev = 10000;
	unsigned int seed = 0;

	bool nonExponential = false;

	double p = 1380.;				// GeV
	double t_min = 0.04;			// GeV^2
	double t_max = 0.7;				// GeV^2
	double th_y_min_cut = 150E-6;	// rad
	double th_y_min_count = 220E-6;	// rad

	double si_th_xy = 121.6E-6;
	double si_de_th_x = 17.43E-6;
	double si_de_th_y = 22.32E-6;

	double de_th_y = 0E-6;

	// parse command line
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
		{
			PrintUsage();
			return 0;
		}

		if (strcmp(argv[i], "--diagonal") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--diagonal' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			diagonal = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "--outputDir") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--outputDir' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			outputDir = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "--de_th_y") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--de_th_y' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			de_th_y = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "--nEv") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--nEv' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			N_ev = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "--seed") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--seed' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			seed = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "--nonExponential") == 0)
		{
			nonExponential = true;
			continue;
		}

		printf("ERROR: unknown option '%s'\n", argv[i]);
		return 1;
	}

	// print configuration
	printf("diagonal = %s\n", diagonal.c_str());
	printf("outputDir = %s\n", outputDir.c_str());
	printf("si_th_xy = %E\n", si_th_xy);
	printf("de_th_y = %E\n", de_th_y);
	printf("N_ev = %u\n", N_ev);
	printf("seed = %u\n", seed);
	printf("nonExponential = %u\n", nonExponential);

	// validate input
	enum { dUnknown, d45b56t, d45t56b } dgn = dUnknown;

	if (diagonal == "45b_56t")
		dgn = d45b56t;

	if (diagonal == "45t_56b")
		dgn = d45t56b;

	if (dgn == dUnknown)
	{
		printf("ERROR: unknown diagonal '%s'.\n", diagonal.c_str());
		return 2;
	}

	// prepare output
	TFile *f_out = TFile::Open((outputDir + "/distributions_" + diagonal + ".root").c_str(), "recreate");

	// prepare plots
	TH1D *h_th_x = new TH1D("h_th_x", ";th_x", 500, -1000E-6, +1000E-6);
	TH1D *h_th_x_L = new TH1D(*h_th_x); h_th_x_L->SetName("h_th_x_L");
	TH1D *h_th_x_R = new TH1D(*h_th_x); h_th_x_R->SetName("h_th_x_R");

	TH1D *h_th_y = new TH1D("h_th_y", ";th_y", 500,
		(dgn == d45t56b) ? -1000E-6 : +100E-6,
		(dgn == d45t56b) ? -100E-6 : +1000E-6
	);
	TH1D *h_th_y_L = new TH1D(*h_th_y); h_th_y_L->SetName("h_th_y_L");
	TH1D *h_th_y_R = new TH1D(*h_th_y); h_th_y_R->SetName("h_th_y_R");

	TGraph *g_th_y_L_F_vs_th_x_L = new TGraph();
	TGraph *g_th_y_L_N_vs_th_x_L = new TGraph();
	TGraph *g_th_y_R_N_vs_th_x_R = new TGraph();
	TGraph *g_th_y_R_F_vs_th_x_R = new TGraph();


	// prepare simulation
	gRandom->SetSeed(seed);

	double la = 1. / 2. / p/p / si_th_xy / si_th_xy;
	double E_min = exp(- la * t_min);
	double E_min_max = exp(- la * t_min) - exp(- la * t_max);
	
	TSpline *icdf_t = NULL;
	if (nonExponential)
		icdf_t = BuildICDF(t_min, t_max);

	// reset counters
	unsigned int n_ev_safe_L_F = 0;
	unsigned int n_ev_safe_L_N = 0;
	unsigned int n_ev_safe_R_N = 0;
	unsigned int n_ev_safe_R_F = 0;

	// run simulation
	unsigned int n_ev_acc = 0;
	while (true)
	{
		// generate event
		double r_t = gRandom->Rndm();
		double t = (nonExponential) ? icdf_t->Eval(r_t) : - log(E_min - r_t * E_min_max) / la;

		double r_phi = gRandom->Rndm();
		double phi = r_phi * M_PI;

		double th = sqrt(t) / p;
		double th_x = th * cos(phi);
		double th_y = th * sin(phi);

		if (th_y < th_y_min_cut)
			continue;

		if (th_y > th_y_min_count)
		  n_ev_acc++;

		if (dgn == d45t56b)
			th_y = -th_y;

		// apply misalignemnt
		th_y += de_th_y;

		// apply smearing
		double sm_th_x_L = si_de_th_x * gRandom->Gaus();
		double sm_th_x_R = si_de_th_x * gRandom->Gaus();

		double th_x_L = th_x + sm_th_x_L;
		double th_x_R = th_x + sm_th_x_R;

		double sm_th_y_L = si_de_th_y * gRandom->Gaus();
		double sm_th_y_R = si_de_th_y * gRandom->Gaus();

		double th_y_L = th_y + sm_th_y_L;
		double th_y_R = th_y + sm_th_y_R;

		double th_y_L_F = th_y_L;
		double th_y_L_N = th_y_L;
		double th_y_R_N = th_y_R;
		double th_y_R_F = th_y_R;

		// update counters
		if (fabs(th_y_L_F) > th_y_min_count) n_ev_safe_L_F++;
		if (fabs(th_y_L_N) > th_y_min_count) n_ev_safe_L_N++;
		if (fabs(th_y_R_N) > th_y_min_count) n_ev_safe_R_N++;
		if (fabs(th_y_R_F) > th_y_min_count) n_ev_safe_R_F++;

		// fill plots
		h_th_x->Fill(th_x);
		h_th_x_L->Fill(th_x_L);
		h_th_x_R->Fill(th_x_R);

		h_th_y->Fill(th_y);
		h_th_y_L->Fill(th_y_L);
		h_th_y_R->Fill(th_y_R);

		int idx = g_th_y_L_F_vs_th_x_L->GetN();
		g_th_y_L_F_vs_th_x_L->SetPoint(idx, th_x_L, th_y_L_F);
		g_th_y_L_N_vs_th_x_L->SetPoint(idx, th_x_L, th_y_L_N);
		g_th_y_R_N_vs_th_x_R->SetPoint(idx, th_x_R, th_y_R_N);
		g_th_y_R_F_vs_th_x_R->SetPoint(idx, th_x_R, th_y_R_F);

		// shall stop ?
		if (n_ev_acc >= N_ev)
			break;
	}

	// print counters
	printf("\n");

	printf("n_ev_safe_L_F = %u\n", n_ev_safe_L_F);
	printf("n_ev_safe_L_N = %u\n", n_ev_safe_L_N);
	printf("n_ev_safe_R_N = %u\n", n_ev_safe_R_N);
	printf("n_ev_safe_R_F = %u\n", n_ev_safe_R_F);

	// save plots
	gDirectory = f_out->mkdir("selected - angles");
	h_th_x->Write();
	h_th_x_L->Write();
	h_th_x_R->Write();

	h_th_y->Write();
	h_th_y_L->Write();
	h_th_y_R->Write();

	g_th_y_L_F_vs_th_x_L->Write("g_th_y_L_F_vs_th_x_L");
	g_th_y_L_N_vs_th_x_L->Write("g_th_y_L_N_vs_th_x_L");
	g_th_y_R_N_vs_th_x_R->Write("g_th_y_R_N_vs_th_x_R");
	g_th_y_R_F_vs_th_x_R->Write("g_th_y_R_F_vs_th_x_R");

	// clean up
	delete f_out;

	return 0;
}
