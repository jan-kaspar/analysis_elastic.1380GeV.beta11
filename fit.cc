#include "../common_definitions.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"

#include <cmath>
#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

vector<TF1 *> ff_list;

//----------------------------------------------------------------------------------------------------

void MakeFit(const string &bt, const string &label, TH1D *h_input)
{
	printf("\n\n================================================== %s ==================================================\n\n",
		(bt + ": " + label).c_str());

	TDirectory *baseDir = gDirectory;
	TDirectory *dataDir = gDirectory = baseDir->mkdir(label.c_str());

	TH1D *h = new TH1D(*h_input);
	h->SetName("dsdt");

	// remove low |t| bins
	for (int bi = 1; bi <= h->GetNbinsX(); ++bi) {
		if (h->GetBinCenter(bi) < anal.t_min_fit) {
			h->SetBinContent(bi, 0.);
			h->SetBinError(bi, 0.);
		}
	}
	
	h->Write();
	
	// make fits for all functions
	for (unsigned int fi = 0; fi < ff_list.size(); fi++) {
		TF1 *ff = new TF1(* ff_list[fi]);

		printf("\n--------- fit %s --------\n", ff->GetName());

		gDirectory = dataDir->mkdir(ff->GetName());

		double t_max_fit = 0.2;
		h->Fit(ff, "Q", "", anal.t_min_fit, t_max_fit);
		h->Fit(ff, "QI", "", anal.t_min_fit, t_max_fit);
		h->SetName("fit");
		h->Write();
	
		// TODO: generalize: a = ff(0), etc.
		double a = ff->GetParameter(0), a_err = ff->GetParError(0);
		double B = ff->GetParameter(1), B_err = ff->GetParError(1);
		double chi2 = ff->GetChisquare(), ndf = ff->GetNDF();
	
		printf("\trange: %.4E to %.4E\n", anal.t_min_fit, t_max_fit);
		printf("\tchi^2 = %.3E, ndf = %.3E, chi^2/ndf = %.3E\n", chi2, ndf, chi2/ndf);
		printf("\ta = %.3E +- %.3E\n", a, a_err);
		printf("\tB = %.3E +- %.3E\n", B, B_err);
		printf("\n");
	
		// integrate cross section
		double t_min_integ = anal.t_min_fit;
		double t_max_integ = 0.4;
		int bi_min = h->FindBin(t_min_integ), bi_max = h->FindBin(t_max_integ);
	
		printf("* summing cross section over bins from %i (left edge %.2E) to %i (right edge %.2E)\n",
			bi_min, h->GetBinLowEdge(bi_min), bi_max, h->GetBinLowEdge(bi_max) + h->GetBinWidth(bi_max));
	
		double integ_sum = 0., integ_sum_err_sq = 0.;
		for (int bi = bi_min; bi <= bi_max; bi++) {
			double w = h->GetBinWidth(bi);
			integ_sum += h->GetBinContent(bi) * w;
			integ_sum_err_sq += pow(h->GetBinError(bi) * w, 2.);
		}
	
		printf("\t%.2f +- %.2f mb\n", integ_sum, sqrt(integ_sum_err_sq));
	
		double t_max_extra = h->GetBinLowEdge(bi_min);
		printf("* integrating between %.2E, %.2E\n", 0., t_max_extra);
		
		/*
		double integ_extra = a/B * (exp(B*t_max_extra) - 1.);
		double cB = a/B*t_max_extra*exp(B*t_max_extra) - integ_extra/B;
		double integ_extra_err = sqrt(integ_extra*integ_extra/a/a*a_err*a_err + cB*cB*B_err*B_err);
		*/

		double integ_extra = ff->Integral(0., t_max_extra);
		double integ_extra_err = 0.; // TODO

		printf("\t%.2f +- %.2f mb\n", integ_extra, integ_extra_err);
	
		double integ_full = integ_sum + integ_extra;
		double integ_full_err = sqrt(integ_sum_err_sq + integ_extra_err*integ_extra_err);
		printf("* full cross section\n\t%.2f +- %.2f mb\n", integ_full, integ_full_err);
	}

	gDirectory = baseDir;
}

//----------------------------------------------------------------------------------------------------

void LoadData(const string &bt, TH1D* &h_t_45b, TH1D* &h_t_45t, TH1D* &h_t_comb)
{
	TFile *fMerged = TFile::Open("merged.root");
	if (fMerged) {
		h_t_45b = (TH1D *) fMerged->Get((bt + "/45 bottom - 56 top/dsdt").c_str());
		h_t_45t = (TH1D *) fMerged->Get((bt + "/45 top - 56 bottom/dsdt").c_str());
		h_t_comb = (TH1D *) fMerged->Get((bt + "/combined/dsdt").c_str());
	} else {
		char hist_name[100];
		sprintf(hist_name, "normalization/h_t_%s_normalized", bt.c_str());
		
		char corr_name[100];
		sprintf(corr_name, "cf,%s/exp3+exp2/corr_final", bt.c_str());

		// get input
		TFile *f_dist_45b = new TFile("distributions_45b_56t.root");
		h_t_45b = (TH1D *) f_dist_45b->Get(hist_name);
		//TFile *f_unf_45b = new TFile("unfolding_45b_56t.root");
		//TH1D *unsm_corr_45b = (TH1D *) f_unf_45b->Get(corr_name);
		//printf("smearing correction for 45b-56t: %p\n", unsm_corr_45b);
		//h_t_45b->Multiply(unsm_corr_45b);
		
		TFile *f_dist_45t = new TFile("distributions_45t_56b.root");
		h_t_45t = (TH1D *) f_dist_45t->Get(hist_name);
		//TFile *f_unf_45t = new TFile("unfolding_45t_56b.root");
		//TH1D *unsm_corr_45t = (TH1D *) f_unf_45t->Get(corr_name);
		//printf("smearing correction for 45t-56b: %p\n", unsm_corr_45t);
		//h_t_45t->Multiply(unsm_corr_45t);
		
		// merge diagonals
		h_t_comb = new TH1D(*h_t_45b);
		h_t_comb->SetName("output");
		for (int i = 1; i <= h_t_comb->GetNbinsX(); i++) {
			double v1 = h_t_45b->GetBinContent(i), e1 = h_t_45b->GetBinError(i);
			double v2 = h_t_45t->GetBinContent(i), e2 = h_t_45t->GetBinError(i);
	
			double w1 = (e1 > 0.) ? 1./e1/e1 : 0.;
			double w2 = (e2 > 0.) ? 1./e2/e2 : 0.;
	
			double v = (w1+w2 > 0.) ? (w1*v1 + w2*v2) / (w1 + w2) : 0.;
			double e = (w1+w2 > 0.) ? 1. / sqrt(w1 + w2) : 0.;
	
			h_t_comb->SetBinContent(i, v);
			h_t_comb->SetBinError(i, e);
	
			//printf("t=%E, v1=%E, v2=%E, v=%E\n", h_t_45b->GetBinCenter(i), v1, v2, v);
		}
	}
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;
	
	Init(argv[1]);
	if (diagonal != dCombined)
		return rcIncompatibleDiagonal;

	// build list of fit functions
	ff_list.push_back(new TF1("exp1", "[0]*exp([1]*x)"));
	ff_list.push_back(new TF1("exp2", "[0]*exp([1]*x + [2]*x*x)"));

	// list of binnings
	vector<string> binnings;
	binnings.push_back("eb");
	binnings.push_back("ub");
	
	// prepare output
	TFile *outF = new TFile("fit.root", "recreate");

	for (unsigned int bi = 0; bi < binnings.size(); bi++) {
		TH1D *h_t_45b, *h_t_45t, *h_t_comb;
		LoadData(binnings[bi], h_t_45b, h_t_45t, h_t_comb);
		
		gDirectory = outF->mkdir(binnings[bi].c_str());

		// make fits
		MakeFit(binnings[bi], "45 bottom - 56 top", h_t_45b);
		MakeFit(binnings[bi], "45 top - 56 bottom", h_t_45t);
		MakeFit(binnings[bi], "combined", h_t_comb);
	}

	delete outF;
	return 0;
}
