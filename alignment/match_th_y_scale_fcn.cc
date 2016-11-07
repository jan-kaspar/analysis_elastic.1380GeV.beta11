#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"

#include <cmath>

#include "kolmogorov.h"

using namespace std;

bool debug;

//----------------------------------------------------------------------------------------------------

double Gauss(double x, double si)
{
	return exp(-x*x / 2. / si/si) / sqrt(2. * M_PI) / si;
}

//----------------------------------------------------------------------------------------------------

double GaussIntegral(double x_min, double x_max, double si)
{
	const double d = si * sqrt(2.);
	return ( TMath::Erf(x_max / d) - TMath::Erf(x_min / d) ) / 2.;
}

//----------------------------------------------------------------------------------------------------

void DoMatch(const string &label, TGraph *g, double th_y_min, double th_y_max, double si_th_ref, FILE *f_res_out, const string &res_label)
{
	printf(">> DoMatch: %s\n", label.c_str());

	// book new subdir
	TDirectory *d_top = gDirectory;
	gDirectory = d_top->mkdir(label.c_str());

	// shift range (mu rad)
	double sh_min = -60E-6, sh_max = +60E-6, sh_step = 1E-6;

	if (debug)
	{
		sh_min = -40E-6;
		sh_max = +40E-6;
		sh_step = 5E-6;
	}

	// book plots
	TGraph *g_kol_meas = new TGraph(); g_kol_meas->SetName("g_kol_meas");

	// shift iterations
	double sh_best_kol = 0;
	double kol_meas_best = 1E100;
	for (double sh = sh_min; sh <= sh_max; sh += sh_step)
	{
		if (debug)
			printf("    sh = %.2f urad\n", sh * 1E6);

		// crop limits
		double th_y_min_sh = th_y_min + sh;
		double th_y_max_sh = th_y_max + sh;

		if (debug)
			printf("    th_y_min_sh = %E, th_y_max_sh = %E\n", th_y_min_sh, th_y_max_sh);

		// make histograms
		unsigned int n_entries = 0;
		CumulativeHistogram ch_th_y;

		TH1D *h_th_y = NULL;
		if (debug)
			h_th_y = new TH1D("", ";th_y", 100, th_y_min_sh, th_y_max_sh);

		for (int i = 0; i < g->GetN(); i++)
		{
			double th_x, th_y;
			g->GetPoint(i, th_x, th_y);

			double th_y_abs = fabs(th_y) + sh;

			if (th_y_min_sh <= th_y_abs && th_y_abs <= th_y_max_sh)
			{
				n_entries++;

				ch_th_y.Fill(th_y_abs);

				if (debug)
					h_th_y->Fill(th_y_abs);
			}
		}

		// run Kolmogorov test
		TGraph *g_kol_diff_graph = NULL, *g_ch_th_x = NULL;
		if (debug)
		{
			g_kol_diff_graph = new TGraph();
			g_ch_th_x = new TGraph();
		}

		double norm = double(n_entries) / GaussIntegral(th_y_min_sh, th_y_max_sh, si_th_ref);
		double sum = 0.;
		double max_diff = -1.;
		for (const auto &p : ch_th_y.data)
		{
			const double &th_y = p.first;
			sum += p.second;

			double F = GaussIntegral(th_y_min_sh, th_y, si_th_ref) * norm;

			double diff = F - sum;

			max_diff = max(max_diff, fabs(diff));

			if (debug)
			{
				g_kol_diff_graph->SetPoint(g_kol_diff_graph->GetN(), th_y, diff);
				g_ch_th_x->SetPoint(g_ch_th_x->GetN(), th_y, F);
			}
		}
		
		double kol_meas = max_diff;

		// fill plots
		int idx = g_kol_meas->GetN();
		g_kol_meas->SetPoint(idx, sh, kol_meas);		

		// debug
		if (debug)
		{
			TDirectory *d_save = gDirectory;
			char buf[100];
			sprintf(buf, "sh = %.2f", sh * 1E6);
			gDirectory = d_save->mkdir(buf);

			TGraph *g_h_th_x = new TGraph();
			g_h_th_x->SetName("g_h_th_x");
			g_h_th_x->SetLineColor(4);
			//for (double th_x = th_y_min_sh; th_x <= th_y_max_sh; th_x += 1E-6)
			for (double th_x = -400E-6; th_x <= +400E-6; th_x += 2E-6)
				g_h_th_x->SetPoint(g_h_th_x->GetN(), th_x, Gauss(th_x, si_th_ref) * norm);

			h_th_y->SetLineColor(2);
			h_th_y->SetName("h_th_y");
			h_th_y->Scale(1., "width");

			TCanvas *c = new TCanvas();
			g_h_th_x->Draw("al");
			h_th_y->Draw("same");
			c->Write("h_cmp");
			delete c;

			g_ch_th_x->SetLineColor(4);
			g_ch_th_x->SetName("g_ch_th_x");
	
			TGraph *g_ch_th_y = ch_th_y.GetGraph();
			g_ch_th_y->SetLineColor(2);
			g_ch_th_y->SetName("g_ch_th_y");

			c = new TCanvas();
			g_ch_th_x->Draw("al");
			g_ch_th_y->Draw("l");
			c->Write("ch_cmp");
			delete c;

			g_kol_diff_graph->SetLineColor(8);
			g_kol_diff_graph->Write("g_kol_diff_graph");
	
			gDirectory = d_save;
		}

		// determine best
		if (kol_meas < kol_meas_best)
		{
			kol_meas_best = kol_meas;
			sh_best_kol = sh;
		}
	}
	// save results
	fprintf(f_res_out, "%s/sh_best_kol,%E\n", res_label.c_str(), sh_best_kol);

	// save plots
	g_kol_meas->SetMarkerStyle(2);
	g_kol_meas->Write();

	// restore global directory
	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: match_th_y_scale_fcn <option> <option> ...\n");
	printf("OPTIONS:\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string inputDir = "../DS1_rel_al_only/";
	string outputDir = "./";
	string outputTag = "match_th_y_scale_fcn";

	double si_th_ref = 123.69E-6; // si_th_x corrected for the different smearing in x and y

	debug = false;

	// parse command line
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "--inputDir") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--inputDir' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			inputDir = argv[++i];
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

		if (strcmp(argv[i], "--outputTag") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--outputTag' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			outputTag = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "--si_th_ref") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--si_th_ref' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			si_th_ref = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "--debug") == 0)
		{
			debug = true;
			continue;
		}

		printf("ERROR: unknown option '%s'\n", argv[i]);
		return 1;
	}

	// print settings
	printf("* inputDir = %s\n", inputDir.c_str());
	printf("* outputDir = %s\n", outputDir.c_str());
	printf("* outputTag = %s\n", outputTag.c_str());
	printf("* si_th_ref = %E\n", si_th_ref);

	// prepare output
	TFile *f_out = TFile::Open((outputDir + "/" + outputTag + ".root").c_str(), "recreate");

	FILE *f_res_out = fopen((outputDir + "/" + outputTag + ".out").c_str(), "w");

	// process input
	//for (const auto &dataset : datasets)
	{
		//TDirectory *d_dataset = f_out->mkdir(dataset.c_str());
		TDirectory *d_dataset = f_out;

		TFile *f_in;
		string diagonal;

		//--------------------

		diagonal = "45b_56t";	
		f_in = TFile::Open((inputDir + "/distributions_" + diagonal + ".root").c_str());
		gDirectory = d_dataset->mkdir(diagonal.c_str());

		DoMatch("L_F", (TGraph *)f_in->Get("selected - angles/g_th_y_L_F_vs_th_x_L"), +220E-6, +420E-6, si_th_ref, f_res_out, "45b_56t/L_F");
		DoMatch("L_N", (TGraph *)f_in->Get("selected - angles/g_th_y_L_N_vs_th_x_L"), +220E-6, +420E-6, si_th_ref, f_res_out, "45b_56t/L_N");
		DoMatch("R_N", (TGraph *)f_in->Get("selected - angles/g_th_y_R_N_vs_th_x_R"), +220E-6, +420E-6, si_th_ref, f_res_out, "45b_56t/R_N");
		DoMatch("R_F", (TGraph *)f_in->Get("selected - angles/g_th_y_R_F_vs_th_x_R"), +220E-6, +420E-6, si_th_ref, f_res_out, "45b_56t/R_F");

		delete f_in;

		//--------------------

		diagonal = "45t_56b";	
		f_in = TFile::Open((inputDir + "/distributions_" + diagonal + ".root").c_str());
		gDirectory = d_dataset->mkdir(diagonal.c_str());

		DoMatch("L_F", (TGraph *)f_in->Get("selected - angles/g_th_y_L_F_vs_th_x_L"), +220E-6, +420E-6, si_th_ref, f_res_out, "45t_56b/L_F");
		DoMatch("L_N", (TGraph *)f_in->Get("selected - angles/g_th_y_L_N_vs_th_x_L"), +220E-6, +420E-6, si_th_ref, f_res_out, "45t_56b/L_N");
		DoMatch("R_N", (TGraph *)f_in->Get("selected - angles/g_th_y_R_N_vs_th_x_R"), +220E-6, +420E-6, si_th_ref, f_res_out, "45t_56b/R_N");
		DoMatch("R_F", (TGraph *)f_in->Get("selected - angles/g_th_y_R_F_vs_th_x_R"), +220E-6, +420E-6, si_th_ref, f_res_out, "45t_56b/R_F");

		delete f_in;
	}

	// clean up
	fclose(f_res_out);

	delete f_out;

	return 0;
}
