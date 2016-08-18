#include "../common_definitions.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"

#include <vector>

using namespace std;

struct Point
{
	double x, y, v, e;
};

vector<Point> points;

//----------------------------------------------------------------------------------------------------

class S2_FCN : public ROOT::Minuit2::FCNBase
{
	public:
		S2_FCN() {}

  		double operator() (const std::vector<double> &) const;
  		double Up() const { return 1.; }

		static double f(double x, double y, const std::vector<double> &);

		TH2D* MakeDiffMap(const TH2D *orig, const std::vector<double> &par);
};

//----------------------------------------------------------------------------------------------------

double S2_FCN::f(double x, double y, const std::vector<double> &par)
{
	const double &A = par[0] * 1E9;
	const double &mu_x = par[1] * 1E-6;
	const double &si_x = par[2] * 1E-6;
	const double &mu_y = par[3] * 1E-6;
	const double &si_y = par[4] * 1E-6;
		
	double de_x = (x - mu_x) / si_x;
	double de_y = (y - mu_y) / si_y;
	double val = A * exp(- (de_x*de_x + de_y*de_y) / 2.);
	
	return val;
}

//----------------------------------------------------------------------------------------------------

double S2_FCN::operator() (const std::vector<double> &par) const
{
	//printf("--------------------------------------------------\n");

	double S2 = 0.;

	for (unsigned int i = 0; i < points.size(); i++) {
		const Point &p = points[i];

		double de_f = (p.v - f(p.x, p.y, par)) / p.e;

		S2 += de_f * de_f;
	}

	//printf("S2 = %E\n", S2);

	return S2;
}

//----------------------------------------------------------------------------------------------------

TH2D* S2_FCN::MakeDiffMap(const TH2D *orig, const std::vector<double> &par)
{
	printf(">> MakeDiffMap\n");

	double bs_x = orig->GetXaxis()->GetBinWidth(1);
	double bs_y = orig->GetYaxis()->GetBinWidth(1);

	double lim_x = 600E-6;	int N_x = ceil(lim_x / bs_x);	lim_x = bs_x * N_x;
	double lim_y = 600E-6;	int N_y = ceil(lim_y / bs_y);	lim_y = bs_y * N_y;

	printf("N_x = %u, lim_x = %E\n", N_x, lim_x);
	printf("N_y = %u, lim_y = %E\n", N_y, lim_y);

	TH2D *h = new TH2D("diff map", "", 2*N_x, -lim_x, +lim_x, 2*N_y, -lim_y, +lim_y);

	for (unsigned int i = 0; i < points.size(); i++) {
		const Point &p = points[i];

		int bi_x = h->GetXaxis()->FindBin(p.x);
		int bi_y = h->GetXaxis()->FindBin(p.y);

		double de_f = (p.v - f(p.x, p.y, par)) / p.e;

		h->SetBinContent(bi_x, bi_y, de_f*de_f);
	}

	return h;
}

//----------------------------------------------------------------------------------------------------

void MakeFit(TH2D *h_45b, TH2D *h_45t)
{
	printf(">> MakeFit\n");
	printf("%p, %p\n", h_45b, h_45t);

	// select points
	vector<TH2D *> hists;
	hists.push_back(h_45b);
	hists.push_back(h_45t);
	for (unsigned int hi = 0; hi < hists.size(); hi++) {
		TH2D *h = hists[hi];
		for (int xi = 1; xi <= h->GetNbinsX(); xi++) {
			for (int yi = 1; yi <= h->GetNbinsY(); yi++) {
				Point p;
				p.x = h->GetXaxis()->GetBinCenter(xi);
				p.y = h->GetYaxis()->GetBinCenter(yi);
				p.v = h->GetBinContent(xi, yi);
				p.e = h->GetBinError(xi, yi);

				// exclude empty bins
				if (p.v == 0. && p.e == 0.)
					continue;

				/*
				// remove partially filled bins at acceptance boundaries
				if (p.x < -46E-6 || p.x > 76E-6 || fabs(p.y) > 98E-6 || (p.y < 14E-6 && p.y > -8E-6))
					continue;
				*/

				// theta cut
				double th = sqrt(p.x*p.x + p.y*p.y);
				if (th < 200E-6 || th > 500E-6)
					continue;

				points.push_back(p);

				//printf("th_x = %.0f, th_y = %.0f, v = %.3E, e = %.3E\n", p.x*1E6, p.y*1E6, p.v, p.e);
			}
		}
	}

	printf("number of points: %lu\n", points.size());

	// initialize fitter
	TFitterMinuit *minuit = new TFitterMinuit();
	S2_FCN fcn;
	minuit->SetMinuitFCN(&fcn);

	// set initial parameters
	minuit->SetParameter(0, "const", 4., 0.1, 0., 0.);	// * 1E9
	minuit->SetParameter(1, "mean x (urad)", 0., 1., 0., 0.);	// in urad
	minuit->SetParameter(2, "sigma x (urad)", 100., 1., 0., 0.);
	minuit->SetParameter(3, "mean y (urad)", 0., 1., 0., 0.);
	minuit->SetParameter(4, "sigma y (urad)", 100., 1., 0., 0.);

	// run fit
	minuit->SetPrintLevel(3);
	minuit->CreateMinimizer();
	minuit->Minimize();

	// 
	vector<double> par;
	for (int i = 0; i < minuit->GetNumberTotalParameters(); i++) {
		par.push_back(minuit->GetParameter(i));
		//printf("%E\n", minuit->GetParameter(i));
	}

	double chiSq = fcn(par);
	double ndf = points.size() - minuit->GetNumberTotalParameters();
	printf("chiSq = %E, ndf = %E, chiSq/ndf = %.3f\n", chiSq, ndf, chiSq/ndf);

	TH2D *m = fcn.MakeDiffMap(h_45b, par);
	m->Write();

	// print results
	/*
	printf("res 0: %E\n", minuit->GetParameter(0));
	printf("res 1: %E\n", minuit->GetParameter(1));
	*/
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
	TH2D *h_45b = (TH2D *) inF_45b->Get("normalization/h_th_y_vs_th_x_normalized");
	
	TFile *inF_45t = new TFile("distributions_45t_56b.root");
	TH2D *h_45t = (TH2D *) inF_45t->Get("normalization/h_th_y_vs_th_x_normalized");
	
	// prepare output
	TFile *outF = new TFile("alignment_final.root", "recreate");

	MakeFit(h_45b, h_45t);

	delete outF;
	return 0;
}
