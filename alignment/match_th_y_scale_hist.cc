#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"

#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

double Gauss(double x, double si)
{
	return exp( - x*x / 2. / si/si ) / sqrt(2. * M_PI) / si;
}

//----------------------------------------------------------------------------------------------------

double GaussIntegral(double x_min, double x_max, double si)
{
	const double d = si * sqrt(2.);
	return ( TMath::Erf(x_max / d) - TMath::Erf(x_min / d) ) / 2.;
}

//----------------------------------------------------------------------------------------------------

TF1 *ff_gauss = new TF1("ff_gauss", "[0] * exp(-x*x / 2 / [1] / [1])");

//----------------------------------------------------------------------------------------------------

void DoMatch(const string &label, TGraph *g, double th_y_min, double th_y_max)
{
	printf(">> DoMatch: %s\n", label.c_str());

	// book new subdir
	TDirectory *d_top = gDirectory;
	gDirectory = d_top->mkdir(label.c_str());

	// shift range (mu rad)
	double sh_min = -50E-6, sh_max = +50E-6, sh_step = 20E-6;

	// book plots
	TGraph *g_chi_sq = new TGraph(); g_chi_sq->SetName("g_chi_sq");
	TGraph *g_chi_sq_bins = new TGraph(); g_chi_sq_bins->SetName("g_chi_sq_bins");
	TGraph *g_chi_sq_norm = new TGraph(); g_chi_sq_norm->SetName("g_chi_sq_norm");

	// shift iterations
	double sh_best = 0;
	double S2_best = 1E100;
	for (double sh = sh_min; sh <= sh_max; sh += sh_step)
	{
		printf("    sh = %+.2E\n", sh);

		// crop limits
		double th_y_sh_min = th_y_min + sh;
		double th_y_sh_max = th_y_max + sh;

		printf("        th_y_sh_min = %E, th_y_sh_max = %E\n", th_y_sh_min, th_y_sh_max);

		// get normalisation
		double norm_th_y = GaussIntegral(th_y_sh_min, th_y_sh_max, 124.4E-6);

		printf("        norm_th_y = %E\n", norm_th_y);

		// make histograms
		unsigned int N_bins = 20;
		TH1D *h_th_x = new TH1D("", ";th_x", N_bins, th_y_sh_min, th_y_sh_max); h_th_x->Sumw2();
		TH1D *h_th_y = new TH1D("", ";th_y", N_bins, th_y_sh_min, th_y_sh_max); h_th_y->Sumw2();

		TH1D *h_th_x_full = new TH1D("", ";th_x", 100, -400E-6, +400E-6); h_th_x_full->Sumw2();

		for (int i = 0; i < g->GetN(); i++)
		{
			double th_x, th_y;
			g->GetPoint(i, th_x, th_y);
			th_y += sh;
			
			if (th_y_sh_min <= th_y && th_y <= th_y_sh_max)
			{
				h_th_x->Fill(th_x);
				h_th_y->Fill(th_y, norm_th_y);

				h_th_x_full->Fill(th_x);
			}
		}

		// normalise histograms
		h_th_x->Scale(1., "width");
		h_th_y->Scale(1., "width");
		h_th_x_full->Scale(1., "width");

		ff_gauss->SetParameters(2E7, 125E-6);
		h_th_x_full->Fit(ff_gauss);

		// calculate chi^2
		double S2 = 0.;
		unsigned int S2_bins = 0;
		for (int bi = 1; bi <= h_th_x->GetNbinsX(); bi++)
		{
			const double &v_x = h_th_x->GetBinContent(bi);
			const double &v_x_unc = h_th_x->GetBinError(bi);
			const double &v_y = h_th_y->GetBinContent(bi);
			const double &v_y_unc = h_th_y->GetBinError(bi);

			if (v_x > 5. && v_y > 5.)
			{
				S2_bins++;
				double unc_sq = v_x_unc*v_x_unc + v_y_unc*v_y_unc;
				double diff = v_x - v_y;
				S2 += diff*diff / unc_sq;
			}
		}
		double S2_norm = (S2_bins > 3) ? S2 / S2_bins : 10.;

		// fill plots
		int idx = g_chi_sq->GetN();
		g_chi_sq->SetPoint(idx, sh, S2);
		g_chi_sq_bins->SetPoint(idx, sh, S2_bins);
		g_chi_sq_norm->SetPoint(idx, sh, S2_norm);

		// debug
		if (true)
		{
			TDirectory *d_save = gDirectory;
			char buf[100];
			sprintf(buf, "sh = %.3E", sh);
			gDirectory = d_save->mkdir(buf);

			h_th_x_full->Write("h_th_x_full");

			ff_gauss->SetLineColor(6);
			ff_gauss->Write("ff_gauss");
	
			TCanvas *c = new TCanvas();
			h_th_x->SetName("h_th_x"); h_th_x->SetLineColor(4); h_th_x->SetTitle(";th_x (blue), th_y (red)"); h_th_x->Draw("");
			h_th_y->SetName("h_th_y"); h_th_y->SetLineColor(2); h_th_y->Draw("same");
			ff_gauss->Draw("same");
			c->Write("cmp");
			delete c;
	
			gDirectory = d_save;
		}

		// release histograms
		delete h_th_y;
		delete h_th_x;
	}

	// save plots
	g_chi_sq->Write();
	g_chi_sq_bins->Write();
	g_chi_sq_norm->Write();

	// restore global directory
	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// settings
	vector<string> datasets = {
		"DS1",
		//"DS1_no_add_alignment",
	};
	
	// prepare output
	TFile *f_out = TFile::Open("match_th_y_fixed_scale.root", "recreate");

	// process input
	for (const auto &dataset : datasets)
	{
		TDirectory *d_dataset = f_out->mkdir(dataset.c_str());

		TFile *f_in;
		string diagonal;

		//--------------------

		diagonal = "45b_56t";	
		f_in = TFile::Open(("../" + dataset + "/distributions_" + diagonal + ".root").c_str());
		gDirectory = d_dataset->mkdir(diagonal.c_str());

		DoMatch("L_F", (TGraph *)f_in->Get("selected - angles/g_th_y_L_F_vs_th_x_L"), +220E-6, +420E-6);
		/*
		DoMatch("L_N", (TGraph *)f_in->Get("selected - angles/g_th_y_L_N_vs_th_x_L"), +220E-6, +420E-6);
		DoMatch("R_N", (TGraph *)f_in->Get("selected - angles/g_th_y_R_N_vs_th_x_R"), +220E-6, +420E-6);
		DoMatch("R_F", (TGraph *)f_in->Get("selected - angles/g_th_y_R_F_vs_th_x_R"), +220E-6, +420E-6);
		*/

		delete f_in;

		//--------------------

		diagonal = "45t_56b";	
		f_in = TFile::Open(("../" + dataset + "/distributions_" + diagonal + ".root").c_str());
		gDirectory = d_dataset->mkdir(diagonal.c_str());

		DoMatch("L_F", (TGraph *)f_in->Get("selected - angles/g_th_y_L_F_vs_th_x_L"), -420E-6, -220E-6);
		/*
		DoMatch("L_N", (TGraph *)f_in->Get("selected - angles/g_th_y_L_N_vs_th_x_L"), -420E-6, -220E-6);
		DoMatch("R_N", (TGraph *)f_in->Get("selected - angles/g_th_y_R_N_vs_th_x_R"), -420E-6, -220E-6);
		DoMatch("R_F", (TGraph *)f_in->Get("selected - angles/g_th_y_R_F_vs_th_x_R"), -420E-6, -220E-6);
		*/

		delete f_in;
	}

	// clean up
	delete f_out;

	return 0;
}
