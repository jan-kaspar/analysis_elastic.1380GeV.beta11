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

double GaussCDF(double v, double si)
{
	return (1. + TMath::Erf(v / si / sqrt(2.))) / 2.;
}

//----------------------------------------------------------------------------------------------------

TSpline* BuildGaussICDF(double si, double v_min, double v_max)
{
	// make graph of CDF
	TGraph *g = new TGraph();
	double cdf_min = -1.;
	double cdf_max = -1.;
	for (double v = v_min; ; v += 1E-6)
	{
		double cdf = GaussCDF(v, si);

		g->SetPoint(g->GetN(), cdf, v);

		if (cdf < cdf_min || cdf_min < 0.)
			cdf_min = cdf;

		if (cdf > cdf_max || cdf_max < 0.)
			cdf_max = cdf;

		if (v >= v_max)
			break;
	}

	// normalise CDF
	for (int i = 0; i < g->GetN(); i++)
	{
		double cdf, v;
		g->GetPoint(i, cdf, v);
		g->SetPoint(i, (cdf - cdf_min) / (cdf_max - cdf_min), v);

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

}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string diagonal = "";
	string outputDir = "./";
	double si_th_x = 120E-6;
	double si_th_y = 120E-6;
	double de_th_y = 0E-6;
	unsigned int N_ev = 10000;
	unsigned int seed = 0;

	// parse command line
	for (int i = 1; i < argc; i++)
	{
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

		printf("ERROR: unknown option '%s'\n", argv[i]);
		return 1;
	}

	// print configuration
	printf("diagonal = %s\n", diagonal.c_str());
	printf("outputDir = %s\n", outputDir.c_str());
	printf("si_th_x = %E\n", si_th_x);
	printf("si_th_y = %E\n", si_th_y);
	printf("de_th_y = %E\n", de_th_y);
	printf("N_ev = %u\n", N_ev);
	printf("seed = %u\n", seed);

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
	TH1D *h_th_x = new TH1D("h_th_x", ";th_x", 100, -1000E-6, +1000E-6);
	TH1D *h_th_y = new TH1D("h_th_y", ";th_y", 100,
		(dgn == d45t56b) ? -1000E-6 : +200E-6,
		(dgn == d45t56b) ? -200E-6 : +1000E-6
	);

	TGraph *g_th_y_L_F_vs_th_x_L = new TGraph();
	TGraph *g_th_y_L_N_vs_th_x_L = new TGraph();
	TGraph *g_th_y_R_N_vs_th_x_R = new TGraph();
	TGraph *g_th_y_R_F_vs_th_x_R = new TGraph();

	// prepare inverse cumulative distribution functions
	TSpline *icdf_th_x = BuildGaussICDF(si_th_x, -700E-6, +700E-6);
	TSpline *icdf_th_y = BuildGaussICDF(si_th_y, +200E-6, +700E-6);

	// prepare simulation
	gRandom->SetSeed(seed);

	// run simulation
	for (unsigned int evi = 0; evi < N_ev; evi++)
	{
		// generate event
		double p_th_x = gRandom->Rndm();
		double th_x = icdf_th_x->Eval(p_th_x);

		//printf("p_th_x = %.3f, th_x = %.3E\n", p_th_x, th_x);

		double p_th_y = gRandom->Rndm();
		double th_y = icdf_th_y->Eval(p_th_y);
		if (dgn == d45t56b)
			th_y = -th_y;

		// apply misalignemnt
		th_y += de_th_y;

		// simple simulation
		double th_x_L = th_x;
		double th_x_R = th_x;

		double th_y_L_F = th_y;
		double th_y_L_N = th_y;
		double th_y_R_N = th_y;
		double th_y_R_F = th_y;

		// fill plots
		h_th_x->Fill(th_x);
		h_th_y->Fill(th_y);

		int idx = g_th_y_L_F_vs_th_x_L->GetN();
		g_th_y_L_F_vs_th_x_L->SetPoint(idx, th_x_L, th_y_L_F);
		g_th_y_L_N_vs_th_x_L->SetPoint(idx, th_x_L, th_y_L_N);
		g_th_y_R_N_vs_th_x_R->SetPoint(idx, th_x_R, th_y_R_N);
		g_th_y_R_F_vs_th_x_R->SetPoint(idx, th_x_R, th_y_R_F);
	}

	// save plots
	gDirectory = f_out->mkdir("selected - angles");
	h_th_x->Write();
	h_th_y->Write();

	g_th_y_L_F_vs_th_x_L->Write("g_th_y_L_F_vs_th_x_L");
	g_th_y_L_N_vs_th_x_L->Write("g_th_y_L_N_vs_th_x_L");
	g_th_y_R_N_vs_th_x_R->Write("g_th_y_R_N_vs_th_x_R");
	g_th_y_R_F_vs_th_x_R->Write("g_th_y_R_F_vs_th_x_R");

	// clean up
	delete f_out;

	return 0;
}
