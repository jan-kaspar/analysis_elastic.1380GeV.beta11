#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TVirtualFitter.h"

#include <cmath>

//----------------------------------------------------------------------------------------------------

TGraphErrors* Crop(TH1D *h, double th_y_cut)
{
	TGraphErrors *g = new TGraphErrors();

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

	return g;
}

//----------------------------------------------------------------------------------------------------

unsigned GetOneChiSq_N;

double GetOneChiSq(TGraphErrors *g, double sh, double A, double si)
{
	double S2 = 0.;
	GetOneChiSq_N = 0;

	for (int i = 0; i < g->GetN(); i++)
	{
		double x, y;
		g->GetPoint(i, x, y);
		double y_u = g->GetErrorY(i);

		if (y_u < 1.5)
			continue;

		double x_sh = x + sh;

		double f = A * exp(- x_sh*x_sh/2./si/si);

		double de = (y - f) / y_u;

		S2 += de * de;
		GetOneChiSq_N++;
	}

	return S2;
}

//----------------------------------------------------------------------------------------------------

TGraphErrors *g_45b_L;
TGraphErrors *g_45b_R;
TGraphErrors *g_45t_L;
TGraphErrors *g_45t_R;

void ObjectiveFunction(int &, double *, double &f, double *par, int )
{
	printf(">> ObjectiveFunction\n");

	double A_45b = par[4];
	double A_45t = par[5];
	double si = par[6];

	printf("\tsh_45b_L = %.2E, sh_45b_R = %.2E, sh_45t_L = %.2E, sh_45t_R = %.2E; A_45b = %.2E, A_45t = %.2E, si = %.2E\n",
			par[0], par[1], par[2], par[3], A_45b, A_45t, si);

	double S2 = 0.;
	unsigned int n_points = 0;
	S2 += GetOneChiSq(g_45b_L, par[0], A_45b, si); n_points += GetOneChiSq_N;
	S2 += GetOneChiSq(g_45b_R, par[1], A_45b, si); n_points += GetOneChiSq_N;
	S2 += GetOneChiSq(g_45t_L, par[2], A_45t, si); n_points += GetOneChiSq_N;
	S2 += GetOneChiSq(g_45t_R, par[3], A_45t, si); n_points += GetOneChiSq_N;

	printf("\t=> S2 = %E, points = %u, S2 / points = %E\n", S2, n_points, S2 / n_points);

	f = S2;
}

//----------------------------------------------------------------------------------------------------

void RunOptimisation()
{
	TF1 *fake_fit_fcn = new TF1("ff", "[0]+[1]+[2]+[3]+[4]+[5]+[6]");
	fake_fit_fcn->SetParameters(0., 0., 0., 0., 2E3, 2E3, 110E-6);
	
	fake_fit_fcn->FixParameter(6, 125E-6);

	/*
	fake_fit_fcn->SetParLimits(0, -10E-6, +10E-6);
	fake_fit_fcn->SetParLimits(1, -10E-6, +10E-6);
	fake_fit_fcn->SetParLimits(2, -10E-6, +10E-6);
	fake_fit_fcn->SetParLimits(3, -10E-6, +10E-6);
	*/

	TGraph *fake_fit_graph = new TGraph();
	fake_fit_graph->SetPoint(0, 1., 2.);
	TVirtualFitter::Fitter(fake_fit_graph)->SetFCN(ObjectiveFunction);

	fake_fit_graph->Fit(fake_fit_fcn, "U", "");

	printf("shift in 45b_L: %.1f +- %.1f urad\n", fake_fit_fcn->GetParameter(0)*1E6, fake_fit_fcn->GetParError(0)*1E6);
	printf("shift in 45b_R: %.1f +- %.1f urad\n", fake_fit_fcn->GetParameter(1)*1E6, fake_fit_fcn->GetParError(1)*1E6);
	printf("shift in 45t_L: %.1f +- %.1f urad\n", fake_fit_fcn->GetParameter(2)*1E6, fake_fit_fcn->GetParError(2)*1E6);
	printf("shift in 45t_R: %.1f +- %.1f urad\n", fake_fit_fcn->GetParameter(3)*1E6, fake_fit_fcn->GetParError(3)*1E6);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	TFile *f_in_45b_56t = TFile::Open("../DS1/distributions_45b_56t.root");
	TFile *f_in_45t_56b = TFile::Open("../DS1/distributions_45t_56b.root");

	TH1D *h_45b_L = (TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_L");
	TH1D *h_45b_R = (TH1D *)f_in_45b_56t->Get("selected - angles/h_th_y_R");
	TH1D *h_45t_L = (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_L");
	TH1D *h_45t_R = (TH1D *)f_in_45t_56b->Get("selected - angles/h_th_y_R");

	/*
	g_45b_L = Crop(h_45b_L, 190E-6);
	g_45b_R = Crop(h_45b_R, 190E-6);
	g_45t_L = Crop(h_45t_L, 225E-6);
	g_45t_R = Crop(h_45t_R, 210E-6);
	*/

	g_45b_L = Crop(h_45b_L, 200E-6);
	g_45b_R = Crop(h_45b_R, 200E-6);
	g_45t_L = Crop(h_45t_L, 235E-6);
	g_45t_R = Crop(h_45t_R, 220E-6);

	/*
	g_45b_L = Crop(h_45b_L, 210E-6);
	g_45b_R = Crop(h_45b_R, 210E-6);
	g_45t_L = Crop(h_45t_L, 240E-6);
	g_45t_R = Crop(h_45t_R, 230E-6);
	*/

	RunOptimisation();

	//delete f_out;
	delete f_in_45b_56t;
	delete f_in_45t_56b;
	return 0;
}
