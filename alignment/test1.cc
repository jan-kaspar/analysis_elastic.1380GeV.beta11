#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include <cmath>

//----------------------------------------------------------------------------------------------------

double AnalyzeOne(const char *label, TH1D *h, double th_y_cut)
{
	TGraphErrors *g = new TGraphErrors();
	g->SetName(label);

	for (int bi = 1; bi < h->GetNbinsX(); bi++)
	{
		double x = h->GetBinCenter(bi);
		double v = h->GetBinContent(bi);
		double v_u = h->GetBinError(bi);

		if (fabs(x) > th_y_cut)
		{
			int idx = g->GetN();
			g->SetPoint(idx, x, v);
			g->SetPointError(idx, 0., v_u);
		}
	}

	TF1 *ff = new TF1("ff", "[0] * exp(-(x-[1])^2 / 2 / [2]^2)");
	ff->SetParameters(1E3, 0., 100E-6);
	ff->SetParNames("norm", "mean", "sigma");
	ff->FixParameter(2., 125E-6);

	g->Fit(ff);

	g->Write();

	// TODO
	//return ff->GetParameter(1);
	return 0.;
}

//----------------------------------------------------------------------------------------------------

void FitTwo(const char *label, TH1D *h_P, TH1D *h_N, double cut_P, double cut_N, double sh_P, double sh_N)
{
	printf("\n>>FitTwo: %s\n", label);

	TGraphErrors *g = new TGraphErrors();
	g->SetName(label);

	for (int bi = 1; bi < h_P->GetNbinsX(); bi++)
	{
		double x = h_P->GetBinCenter(bi);
		double v = h_P->GetBinContent(bi);
		double v_u = h_P->GetBinError(bi);

		if (fabs(x) > cut_P)
		{
			int idx = g->GetN();
			g->SetPoint(idx, x + sh_P, v);
			g->SetPointError(idx, 0., v_u);
		}
	}

	for (int bi = 1; bi < h_N->GetNbinsX(); bi++)
	{
		double x = h_N->GetBinCenter(bi);
		double v = h_N->GetBinContent(bi);
		double v_u = h_N->GetBinError(bi);

		if (fabs(x) > cut_N)
		{
			int idx = g->GetN();
			g->SetPoint(idx, x + sh_N, v);
			g->SetPointError(idx, 0., v_u);
		}
	}

	TF1 *ff = new TF1("ff", "[0] * exp(-(x-[1])^2 / 2 / [2]^2)");
	ff->SetParNames("norm", "mean", "sigma");
	ff->SetParameters(1E3, 0., 100E-6);
	//ff->FixParameter(2., 110E-6);

	g->Fit(ff);

	g->Write();

	//return ff->GetParameter(1);
}

//----------------------------------------------------------------------------------------------------

double FitOneThx(const char *label, TH1D *h, double th_x_cut)
{
	printf("\n>>FitOneThx: %s\n", label);

	TGraphErrors *g = new TGraphErrors();
	g->SetName(label);

	for (int bi = 1; bi < h->GetNbinsX(); bi++)
	{
		double x = h->GetBinCenter(bi);
		double v = h->GetBinContent(bi);
		double v_u = h->GetBinError(bi);

		if (fabs(x) < th_x_cut)
		{
			int idx = g->GetN();
			g->SetPoint(idx, x, v);
			g->SetPointError(idx, 0., v_u);
		}
	}

	TF1 *ff = new TF1("ff", "[0] * exp(-(x-[1])^2 / 2 / [2]^2)");
	ff->SetParameters(230., 0., 100E-6);
	ff->SetParNames("norm", "mean", "sigma");

	g->Fit(ff);

	g->Write();

	return 0.;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	TFile *f_in_45b_56t = TFile::Open("../DS1/distributions_45b_56t.root");
	TFile *f_in_45t_56b = TFile::Open("../DS1/distributions_45t_56b.root");

	TFile *f_out = TFile::Open("test1.root", "recreate");

	gDirectory = f_out->mkdir("th_y");

	double sh_45b_L = AnalyzeOne("45b_56t, L", (TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_L"), 190E-6);
	double sh_45b_R = AnalyzeOne("45b_56t, R", (TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_R"), 190E-6);

	double sh_45t_L = AnalyzeOne("45t_56b, L", (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_L"), 225E-6);
	double sh_45t_R = AnalyzeOne("45t_56b, R", (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_R"), 210E-6);

	FitTwo("L",
			(TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_L"), (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_L"),
			190E-6, 225E-6,
			sh_45b_L, sh_45t_L);

	FitTwo("R",
			(TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_R"), (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_R"),
			190E-6, 210E-6,
			sh_45b_R, sh_45t_R);

	// -120E-6, +140E-6

	FitTwo("L_F",
			(TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_L_F"), (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_L_F"),
			220E-6, 220E-6, 0., 0.);

	FitTwo("L_N",
			(TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_L_N"), (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_L_N"),
			220E-6, 220E-6, 0., 0.);

	FitTwo("R_N",
			(TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_R_N"), (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_R_N"),
			220E-6, 220E-6, 0., 0.);

	FitTwo("R_F",
			(TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_R_F"), (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_R_F"),
			220E-6, 220E-6, 0., 0.);

	gDirectory = f_out->mkdir("th_x");
	FitOneThx("45b_56t", (TH1D *)f_in_45b_56t->Get("selected - angles/h_th_x"), 300E-6);
	FitOneThx("45t_56b", (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_x"), 300E-6);

	delete f_out;
	delete f_in_45b_56t;
	delete f_in_45t_56b;

	return 0;
}
