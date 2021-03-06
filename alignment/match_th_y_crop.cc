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

void DoMatch(const string &label, TGraph *g, double th_y_min, double th_y_max, FILE *f_res_out, const string &res_label)
{
	printf(">> DoMatch: %s\n", label.c_str());

	// book new subdir
	TDirectory *d_top = gDirectory;
	gDirectory = d_top->mkdir(label.c_str());

	// shift range (mu rad)
	double sh_min = -60E-6, sh_max = +60E-6, sh_step = 1E-6;

	if (debug)
		sh_step = 10E-6;

	// book plots
	TGraph *g_chi_sq = new TGraph(); g_chi_sq->SetName("g_chi_sq");
	TGraph *g_chi_sq_norm = new TGraph(); g_chi_sq_norm->SetName("g_chi_sq_norm");

	TGraph *g_kol_diff_asc = new TGraph(); g_kol_diff_asc->SetName("g_kol_diff_asc");
	TGraph *g_kol_diff_des = new TGraph(); g_kol_diff_des->SetName("g_kol_diff_des");
	TGraph *g_kol_meas_asc = new TGraph(); g_kol_meas_asc->SetName("g_kol_meas_asc");

	// shift iterations
	double sh_best_kol = 0;
	double kol_meas_best = 1E100;
	for (double sh = sh_min; sh <= sh_max; sh += sh_step)
	{
		if (debug)
			printf("    sh = %.2f urad\n", sh * 1E6);

		// crop limits
		double th_xy_min = th_y_min + sh;
		double th_xy_max = th_y_max + sh;

		if (debug)
			printf("    th_xy_min = %E, th_xy_max = %E\n", th_xy_min, th_xy_max);

		// make histograms
		unsigned int N_bins = 20;
		TH1D *h_th_x = new TH1D("", ";th_x", N_bins, th_xy_min, th_xy_max);
		TH1D *h_th_y = new TH1D("", ";th_y", N_bins, th_xy_min, th_xy_max);
	
		CumulativeHistogram ch_th_y, ch_th_x;

		for (int i = 0; i < g->GetN(); i++)
		{
			double th_x, th_y;
			g->GetPoint(i, th_x, th_y);

			double th_x_abs = fabs(th_x);
			double th_y_abs = fabs(th_y) + sh;

			if (th_xy_min <= th_x_abs && th_x_abs <= th_xy_max && th_xy_min <= th_y_abs && th_y_abs <= th_xy_max)
			{
				h_th_x->Fill(th_x_abs);
				h_th_y->Fill(th_y_abs);

				ch_th_x.Fill(th_x_abs);
				ch_th_y.Fill(th_y_abs);
			}
		}

		// calculate chi^2
		double S2 = 0.;
		unsigned int S2_bins = 0;
		for (int bi = 1; bi <= h_th_x->GetNbinsX(); bi++)
		{
			const double &v_x = h_th_x->GetBinContent(bi);
			const double &v_x_unc = h_th_x->GetBinError(bi);
			const double &v_y = h_th_y->GetBinContent(bi);
			const double &v_y_unc = h_th_y->GetBinError(bi);

			if (v_x > 10. && v_y > 10.)
			{
				S2_bins++;
				double unc_sq = v_x_unc*v_x_unc + v_y_unc*v_y_unc;
				double diff = v_x - v_y;
				S2 += diff*diff / unc_sq;
			}
		}
		double S2_norm = (S2_bins > 3) ? S2 / S2_bins : 100.;

		// run Kolmogorov test
		TGraph *g_kol_graph_asc = NULL;
		if (debug)
			g_kol_graph_asc = new TGraph();

		double kol_diff_asc = KolmogorovTest(ch_th_x, ch_th_y, true, g_kol_graph_asc);
		double kol_diff_des = KolmogorovTest(ch_th_x, ch_th_y, false);

		double kol_meas_asc = fabs(kol_diff_asc) / sqrt(h_th_x->GetEntries());

		// fill plots
		int idx = g_chi_sq->GetN();
		g_chi_sq->SetPoint(idx, sh, S2);
		g_chi_sq_norm->SetPoint(idx, sh, S2_norm);
		g_kol_diff_asc->SetPoint(idx, sh, kol_diff_asc);
		g_kol_diff_des->SetPoint(idx, sh, kol_diff_des);
		g_kol_meas_asc->SetPoint(idx, sh, kol_meas_asc);

		// debug
		if (debug)
		{
			TDirectory *d_save = gDirectory;
			char buf[100];
			sprintf(buf, "sh = %.2f", sh * 1E6);
			gDirectory = d_save->mkdir(buf);
	
			TCanvas *c = new TCanvas();
			h_th_x->SetName("h_th_x"); h_th_x->SetLineColor(2); h_th_x->Draw("");
			h_th_y->SetName("h_th_y"); h_th_y->SetLineColor(4); h_th_y->Draw("same");
			c->Write("h_cmp");
			delete c;

			TGraph *g_ch_th_x = ch_th_x.GetGraph();
			g_ch_th_x->SetLineColor(2);
			g_ch_th_x->SetName("g_ch_th_x");

			TGraph *g_ch_th_y = ch_th_y.GetGraph();
			g_ch_th_y->SetLineColor(4);
			g_ch_th_y->SetName("g_ch_th_y");

			c = new TCanvas();
			g_ch_th_x->Draw("al");
			g_ch_th_y->Draw("l");
			c->Write("ch_cmp");
			delete c;

			g_kol_graph_asc->Write("g_kol_graph_asc");
	
			gDirectory = d_save;
		}

		// determine best
		if (kol_meas_asc < kol_meas_best)
		{
			kol_meas_best = kol_meas_asc;
			sh_best_kol = sh;
		}
	
		// release histograms
		delete h_th_y;
		delete h_th_x;
	}

	// print results
	printf("    sh_best_kol = %E\n", sh_best_kol);

	// save results
	fprintf(f_res_out, "%s/sh_best_kol,%E\n", res_label.c_str(), sh_best_kol);

	// save plots
	g_chi_sq->Write();
	g_chi_sq_norm->Write();
	g_kol_diff_asc->Write();
	g_kol_diff_des->Write();
	g_kol_meas_asc->Write();

	// restore global directory
	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: match_th_y_crop <option> <option> ...\n");
	printf("OPTIONS:\n");

}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string inputDir = "../DS1_rel_al_only/";
	string outputDir = "./";
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

		if (strcmp(argv[i], "--debug") == 0)
		{
			debug = true;
			continue;
		}

		printf("ERROR: unknown option '%s'\n", argv[i]);
		return 1;
	}

	// prepare output
	TFile *f_out = TFile::Open((outputDir + "/match_th_y_crop.root").c_str(), "recreate");

	FILE *f_res_out = fopen((outputDir + "/match_th_y_crop.out").c_str(), "w");

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

		DoMatch("L_F", (TGraph *)f_in->Get("selected - angles/g_th_y_L_F_vs_th_x_L"), +220E-6, +420E-6, f_res_out, "45b_56t/L_F");
		// TODO: uncomment
		/*
		DoMatch("L_N", (TGraph *)f_in->Get("selected - angles/g_th_y_L_N_vs_th_x_L"), +220E-6, +420E-6, f_res_out, "45b_56t/L_N");
		DoMatch("R_N", (TGraph *)f_in->Get("selected - angles/g_th_y_R_N_vs_th_x_R"), +220E-6, +420E-6, f_res_out, "45b_56t/R_N");
		DoMatch("R_F", (TGraph *)f_in->Get("selected - angles/g_th_y_R_F_vs_th_x_R"), +220E-6, +420E-6, f_res_out, "45b_56t/R_F");
		*/

		delete f_in;

		//--------------------

		diagonal = "45t_56b";	
		f_in = TFile::Open((inputDir + "/distributions_" + diagonal + ".root").c_str());
		gDirectory = d_dataset->mkdir(diagonal.c_str());

		DoMatch("L_F", (TGraph *)f_in->Get("selected - angles/g_th_y_L_F_vs_th_x_L"), +220E-6, +420E-6, f_res_out, "45t_56b/L_F");
		// TODO: uncomment
		/*
		DoMatch("L_N", (TGraph *)f_in->Get("selected - angles/g_th_y_L_N_vs_th_x_L"), +220E-6, +420E-6, f_res_out, "45t_56b/L_N");
		DoMatch("R_N", (TGraph *)f_in->Get("selected - angles/g_th_y_R_N_vs_th_x_R"), +220E-6, +420E-6, f_res_out, "45t_56b/R_N");
		DoMatch("R_F", (TGraph *)f_in->Get("selected - angles/g_th_y_R_F_vs_th_x_R"), +220E-6, +420E-6, f_res_out, "45t_56b/R_F");
		*/

		delete f_in;
	}

	// clean up
	fclose(f_res_out);

	delete f_out;

	return 0;
}
