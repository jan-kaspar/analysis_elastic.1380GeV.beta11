#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TH1D.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TGraphErrors.h"

#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

bool debug = false;

double th_min = 220E-6;
double th_max = 400E-6;

//----------------------------------------------------------------------------------------------------

class S2_FCN : public ROOT::Minuit2::FCNBase
{
	public:
		S2_FCN() {}

  		double operator() (const std::vector<double> &) const;
  		double Up() const { return 1.; }

		static double f(double x, const std::vector<double> &);
};


//----------------------------------------------------------------------------------------------------

double S2_FCN::f(double x, const std::vector<double> &par)
{
	const double &A = par[0];
	const double &mu = par[1];
	const double &si = par[2];
		
	double de_x = (x - mu) / si;
	double val = A * exp(- de_x*de_x / 2.);
	
	return val;
}

//----------------------------------------------------------------------------------------------------

TH1D *hist_to_fit;
unsigned int points;
bool cut_outliers;

int verbosity = 0;

bool makeRelDiffPlot = false;
TGraphErrors *relDiffPlot = NULL;

double S2_FCN::operator() (const std::vector<double> &par) const
{
	//printf("--------------------------------------------------\n");

	if (makeRelDiffPlot)
		relDiffPlot = new TGraphErrors();

	double S2 = 0.;
	points = 0;
	for (int bi = 1; bi <= hist_to_fit->GetNbinsX(); bi++)
	{
		double c = hist_to_fit->GetBinCenter(bi);
		double v = hist_to_fit->GetBinContent(bi);
		double v_u = hist_to_fit->GetBinError(bi);

		// skip empty bins
		if (fabs(v) < 1E-3)
			continue;

		double de_f = (v - f(c, par)) / v_u;

		if (cut_outliers && fabs(de_f) > 4)
			continue;

		points++;
		S2 += de_f * de_f;

		if (verbosity)
			printf("bi=%i, de_f = %.2f, S2_inc=%E\n", bi, de_f, de_f*de_f);

		if (makeRelDiffPlot)
			relDiffPlot->SetPoint(relDiffPlot->GetN(), c, de_f);
	}

	//printf("S2 = %E\n", S2);

	return S2;
}

//----------------------------------------------------------------------------------------------------

void DoFit(TH1D *h, double &si, double &si_u)
{
	printf("\n>> DoFit(%s)\n", h->GetName());
	hist_to_fit = h;

	// histogram maximum
	double h_max = 0.;
	for (int bi = 1; bi <= h->GetNbinsX(); ++bi)
	{
		double v = h->GetBinContent(bi);
		h_max = max(h_max, v);
	}

	printf("\th_max = %.2E\n", h_max);

	// initialize fitter
	TFitterMinuit *minuit = new TFitterMinuit();
	S2_FCN fcn;
	minuit->SetMinuitFCN(&fcn);

	// set initial parameters
	minuit->SetParameter(0, "const", 1.5*h_max, 0.1*h_max, 0., 0.);
	minuit->SetParameter(1, "mean", 0., 1., 0., 0.);
	minuit->SetParameter(2, "sigma", 40., 1., 0., 0.);

	printf("* angles in urad\n");

	// run fit
	verbosity = 0.;
	makeRelDiffPlot = false;
	minuit->SetPrintLevel(0);
	minuit->CreateMinimizer();

	cut_outliers = false;
	minuit->Minimize();
	
	cut_outliers = true;
	minuit->Minimize();

	// print results
	printf("\tconst = %.3f +- %.3f\n", minuit->GetParameter(0), minuit->GetParError(0));
	printf("\tmean = %.3f +- %.3f\n", minuit->GetParameter(1), minuit->GetParError(1));
	printf("\tsigma = %.3f +- %.3f\n", minuit->GetParameter(2), minuit->GetParError(2));
	
	si = minuit->GetParameter(2);
	si_u = minuit->GetParError(2);

	// get minuit parameters
	vector<double> par;
	for (int i = 0; i < minuit->GetNumberTotalParameters(); i++)
	{
		par.push_back(minuit->GetParameter(i));
	}

	verbosity = 0;
	makeRelDiffPlot = true;
	double chiSq = fcn(par);
	double ndf = points - minuit->GetNumberTotalParameters();
	printf("\tchiSq = %E, ndf = %E, chiSq/ndf = %.3f\n", chiSq, ndf, chiSq/ndf);

	// save rel. diff. plot
	relDiffPlot->SetName((string(h->GetName()) + "_fit_rel_diff").c_str());
	relDiffPlot->SetMarkerStyle(20);
	relDiffPlot->SetMarkerSize(0.7);
	relDiffPlot->Write();
}

//----------------------------------------------------------------------------------------------------


void DoScalingFitTest(TGraph *g_th_45b, TGraph *g_corr_45b=NULL, TGraph *g_th_45t=NULL, TGraph *g_corr_45t=NULL,
		bool th_abs=true)
{
	vector<TGraph *> g_th;
	vector<TGraph *> g_corr;

	g_th.push_back(g_th_45b); g_corr.push_back(g_corr_45b);
	if (g_th_45t)
	{
		g_th.push_back(g_th_45t);
		g_corr.push_back(g_corr_45t);
	}

	TH1D *h_th_x = new TH1D("h_th_x", ";#theta_x   (#murad)", 200, -500., 500.); h_th_x->SetLineColor(2);
	TH1D *h_th_y = new TH1D("h_th_y", ";#theta_y   (#murad)", 200, -500., 500.); h_th_y->SetLineColor(4);

	h_th_x->Sumw2();
	h_th_y->Sumw2();

	TDirectory *topDir = gDirectory;

	// build histograms
	for (unsigned int gi = 0; gi < g_th.size(); gi++)
	{
		double *X = g_th[gi]->GetX(), *Y = g_th[gi]->GetY();
		double *div_corr = NULL, *norm_corr = NULL;
		if (g_corr[gi])
		{
			div_corr = g_corr[gi]->GetX();
			norm_corr = g_corr[gi]->GetY();
		}

		for (int pi = 0; pi < g_th[gi]->GetN(); pi++)
		{
			double x = X[pi], y = Y[pi];

			if (fabs(y) < th_min || fabs(y) > th_max)
				continue;
			
			if (fabs(x) < th_min || fabs(x) > th_max)
				continue;
	
			double w = 1.;
			if (div_corr)
				w = div_corr[pi] * norm_corr[pi];

			if (th_abs)
			{
				h_th_x->Fill(fabs(x)*1E6, w);
				h_th_y->Fill(fabs(y)*1E6, w);
			} else {
				h_th_x->Fill(x*1E6, w);
				h_th_y->Fill(y*1E6, w);
			}
		}
	}

	// fit histograms
	double si_th_x=0., si_th_x_unc=0.;
	DoFit(h_th_x, si_th_x, si_th_x_unc);

	double si_th_y=0., si_th_y_unc=0.;
	DoFit(h_th_y, si_th_y, si_th_y_unc);

	double c = si_th_y / si_th_x;
	double c_unc = sqrt( si_th_y_unc*si_th_y_unc / si_th_x/si_th_x + si_th_y*si_th_y/si_th_x/si_th_x/si_th_x/si_th_x * si_th_x_unc*si_th_x_unc );

	printf("th_x scaling: c = %.3f +- %.3f\n", c, c_unc);
	
	// write histograms
	h_th_x->Write();
	h_th_y->Write();

	// write result
	TGraph *g_c = new TGraph();
	g_c->SetName("g_c");
	g_c->SetPoint(0, c, c_unc);
	g_c->Write();
	
	gDirectory = topDir;

	delete h_th_x;
	delete h_th_y;
}

//----------------------------------------------------------------------------------------------------

double pearson_min = 220.;
double pearson_max = 400.;

double PearsonTest(TH1D *h1, TH1D *h2)
{
	double S = 0.;
	unsigned int N_points = 0;

	for (int bi = 1; bi <= h1->GetNbinsX(); ++bi)
	{
		double c1 = h1->GetBinCenter(bi);
		double v1 = h1->GetBinContent(bi);
		double v_u1 = h1->GetBinError(bi);

		double c2 = h2->GetBinCenter(bi);
		double v2 = h2->GetBinContent(bi);
		double v_u2 = h2->GetBinError(bi);

		// different binning?
		if (fabs(c1 - c2) / (c1 + c2) > 1E-3)
			printf("ERROR: c1 != c2.\n");

		// safe range
		if (fabs(c1) < pearson_min || fabs(c1) > pearson_max)
			continue;

		// skip empty bins
		if (fabs(v1) < 1E-3 && fabs(v2) < 1E-3)
			continue;

		double v_u_eff = sqrt(v_u1*v_u1 + v_u2*v_u2);
		double d = (v1 - v2) / v_u_eff;

		S += d*d;
		N_points++;
	}

	printf("\tS2 = %.2E, N_points = %u, S2 / N_points = %.2f\n", S, N_points, S / N_points);

	return S/N_points;
	//return S;
}

//----------------------------------------------------------------------------------------------------

void DoScalingPearsonTest(TGraph *g_th_45b, TGraph *g_corr_45b=NULL, TGraph *g_th_45t=NULL, TGraph *g_corr_45t=NULL,
		bool th_abs=true)
{
	vector<TGraph *> g_th;
	vector<TGraph *> g_corr;

	g_th.push_back(g_th_45b); g_corr.push_back(g_corr_45b);
	if (g_th_45t)
		{ g_th.push_back(g_th_45t); g_corr.push_back(g_corr_45t); }

	TGraph *g = new TGraph();
	g->SetName("pearson test");
	g->SetTitle(";factor scaling #theta_{x}");

	TH1D *h_th_x = new TH1D("h_th_x", ";#theta_x   (#murad)", 200, -500., 500.); h_th_x->SetLineColor(2);
	TH1D *h_th_y = new TH1D("h_th_y", ";#theta_y   (#murad)", 200, -500., 500.); h_th_y->SetLineColor(4);

	h_th_x->Sumw2();
	h_th_y->Sumw2();

	TDirectory *topDir = gDirectory;

	for (double scale_x = 0.85; scale_x <= 1.35; scale_x += 0.002)
	{
		printf("scale_x = %.4f\n", scale_x);
		
		bool saveDetails = debug || fabs(scale_x - 1.) < 1E-6;

		if (saveDetails)
		{
			char buf[100];
			sprintf(buf, "%.3f", scale_x);
			gDirectory = topDir->mkdir(buf);
		}

		// build histograms
		h_th_x->Reset();
		h_th_y->Reset();

		for (unsigned int gi = 0; gi < g_th.size(); gi++)
		{
			double *X = g_th[gi]->GetX(), *Y = g_th[gi]->GetY();
			double *div_corr = NULL, *norm_corr = NULL;
			if (g_corr[gi])
				{ div_corr = g_corr[gi]->GetX(); norm_corr = g_corr[gi]->GetY(); }

			for (int pi = 0; pi < g_th[gi]->GetN(); pi++)
			{
				double x = X[pi] * scale_x, y = Y[pi];

				if (fabs(y) < th_min || fabs(y) > th_max)
					continue;
				
				if (fabs(x) < th_min || fabs(x) > th_max)
					continue;
		
				double w = 1.;
				if (div_corr)
					w = div_corr[pi] * norm_corr[pi];
		
				if (th_abs)
				{
					h_th_x->Fill(fabs(x)*1E6, w);
					h_th_y->Fill(fabs(y)*1E6, w);
				} else {
					h_th_x->Fill(x*1E6, w);
					h_th_y->Fill(y*1E6, w);
				}
			}
		}
		
		if (saveDetails)
		{
			h_th_x->Write();
			h_th_y->Write();
		}

		// apply test
		double test_value = PearsonTest(h_th_x, h_th_y);

		// add point to graph
		g->SetPoint(g->GetN(), scale_x, test_value);
	}
	
	gDirectory = topDir;

	g->Write();

	delete h_th_x;
	delete h_th_y;
}

//----------------------------------------------------------------------------------------------------

void DoScalingTests(TGraph *g_th_45b, TGraph *g_corr_45b=NULL, TGraph *g_th_45t=NULL, TGraph *g_corr_45t=NULL,
		bool th_abs=true)
{
	TDirectory *topDir = gDirectory;

	gDirectory = topDir->mkdir("fit");
	DoScalingFitTest(g_th_45b, g_corr_45b, g_th_45t, g_corr_45t, th_abs);

	gDirectory = topDir->mkdir("pearson");
	DoScalingPearsonTest(g_th_45b, g_corr_45b, g_th_45t, g_corr_45t, th_abs);
	
	gDirectory = topDir;
}

//----------------------------------------------------------------------------------------------------

void RunGlobalTest(TFile *f_45b, TFile *f_45t)
{
	TGraph *g_th_45b = (TGraph *) f_45b->Get("acceptance correction/g_th_y_vs_th_x_acc");
	TGraph *g_corr_45b = (TGraph *) f_45b->Get("normalization/g_norm_corr_vs_div_corr");

	TGraph *g_th_45t = (TGraph *) f_45t->Get("acceptance correction/g_th_y_vs_th_x_acc");
	TGraph *g_corr_45t = (TGraph *) f_45t->Get("normalization/g_norm_corr_vs_div_corr");

	TDirectory *topDir = gDirectory;

	printf("\n---------- global ----------\n");
	gDirectory = topDir->mkdir("global");
	DoScalingTests(g_th_45b, g_corr_45b, g_th_45t, g_corr_45t, false);

	gDirectory = topDir;
}

//----------------------------------------------------------------------------------------------------

void RunRPPairTests(TFile *f_45b_56t, TFile *f_45t_56b)
{
	TGraph *g_th_45b = (TGraph *) f_45b_56t->Get("selected - angles/g_th_y_L_vs_th_x_L");
	TGraph *g_th_56t = (TGraph *) f_45b_56t->Get("selected - angles/g_th_y_R_vs_th_x_R");

	TGraph *g_th_45t = (TGraph *) f_45t_56b->Get("selected - angles/g_th_y_L_vs_th_x_L");
	TGraph *g_th_56b = (TGraph *) f_45t_56b->Get("selected - angles/g_th_y_R_vs_th_x_R");

	TDirectory *topDir = gDirectory;

	printf("\n---------- 45 top ----------\n");
	gDirectory = topDir->mkdir("45_top");
	DoScalingTests(g_th_45t);
	
	printf("\n---------- 45 bottom ----------\n");
	gDirectory = topDir->mkdir("45_bottom");
	DoScalingTests(g_th_45b);

	printf("\n---------- 56 top ----------\n");
	gDirectory = topDir->mkdir("56_top");
	DoScalingTests(g_th_56t);

	printf("\n---------- 56 bottom ----------\n");
	gDirectory = topDir->mkdir("56_bottom");
	DoScalingTests(g_th_56b);

	gDirectory = topDir;
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal
	Init(argv[1]);
	if (diagonal != dCombined)
		return rcIncompatibleDiagonal;

	// get input data
	TFile *inF_45b = new TFile("distributions_45b_56t.root");
	TFile *inF_45t = new TFile("distributions_45t_56b.root");
	
	// prepare output
	TFile *outF = new TFile("scale_th_x.root", "recreate");

	// run tests
	RunGlobalTest(inF_45b, inF_45t);
	// TODO
	//RunRPPairTests(inF_45b, inF_45t);

	delete outF;
	return 0;
}
