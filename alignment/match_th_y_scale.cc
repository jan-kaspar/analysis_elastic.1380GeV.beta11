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

struct Entry
{
	TH1D *h_th_y;
	double th_y_min, th_y_max;	// rad
	double si_ref;				// rad
};

//----------------------------------------------------------------------------------------------------

double Gauss(double x, double si)
{
	return exp( - x*x / 2. / si/si );
}

//----------------------------------------------------------------------------------------------------

double GaussIntegral(double x_min, double x_max, double si)
{
	const double d = si * sqrt(2.);
	return si * sqrt(M_PI / 2.) * ( TMath::Erf(x_max / d) - TMath::Erf(x_min / d) );
}

//----------------------------------------------------------------------------------------------------

void DoMatch(const vector<Entry> &entries, const string &label = "")
{
	printf("* DoMatch: %s\n", label.c_str());

	TDirectory *d_top = gDirectory;
	gDirectory = d_top->mkdir(label.c_str());

	for (const auto &e : entries)
		printf("    e.h_th_y = %p\n", e.h_th_y);

	// determine auxiliary quantities
	struct EntryAux
	{
		int bin_min, bin_max;
		double integ_h;
	};

	vector<EntryAux> entriesAux;

	for (const auto &e : entries)
	{
		EntryAux aux;

		aux.bin_min = e.h_th_y->GetXaxis()->FindBin(e.th_y_min);
		if (aux.bin_min < 1)
			aux.bin_min = 1;

		aux.bin_max = e.h_th_y->GetXaxis()->FindBin(e.th_y_max);
		if (aux.bin_max >= e.h_th_y->GetNbinsX())
			aux.bin_max = e.h_th_y->GetNbinsX();

		aux.integ_h = 0.;
		for (int bin = aux.bin_min; bin <= aux.bin_max; ++bin)
			aux.integ_h += e.h_th_y->GetBinContent(bin) * e.h_th_y->GetBinWidth(bin);

		entriesAux.push_back(aux);

		printf("    bin_min=%i, bin_max=%i, integ_h=%E\n", aux.bin_min, aux.bin_max, aux.integ_h);
	}

	// book plots
	TGraph *g_chi_sq = new TGraph(); g_chi_sq->SetName("g_chi_sq");
	
	// shift range
	double sh_min = -100E-6, sh_max = +200E-6, sh_step = 10E-6;

	// shift iterations
	double sh_best = 0;
	double S2_best = 1E100;
	vector<double> norm_vec_best;
	for (double sh = sh_min; sh <= sh_max; sh += sh_step)
	{
		printf("    sh = %E\n", sh);

		double S2 = 0.;
		vector<double> norm_vec;
		for (unsigned int ei = 0; ei < entries.size(); ei++)
		{
			const Entry &entry = entries[ei];
			const EntryAux &aux = entriesAux[ei];

			// get normalisation
			double integ_th_y_min = entry.h_th_y->GetBinLowEdge(aux.bin_min);
			double integ_th_y_max = entry.h_th_y->GetBinLowEdge(aux.bin_max) + entry.h_th_y->GetBinWidth(aux.bin_max);

			printf("        integ_th_y_min = %E, integ_th_y_max = %E\n", integ_th_y_min, integ_th_y_max);

			double integ_fcn = GaussIntegral(integ_th_y_min + sh, integ_th_y_max + sh, entry.si_ref);

			printf("        integ_fcn = %E\n", integ_fcn);

			double norm = aux.integ_h / integ_fcn;
			norm_vec.push_back(norm);
		
			// get chi_sq
			for (int bin = aux.bin_min; bin <= aux.bin_max; ++bin)
			{
				const double &x_cen = entry.h_th_y->GetBinCenter(bin) + sh;
				const double &y = entry.h_th_y->GetBinContent(bin);
				const double &y_unc = entry.h_th_y->GetBinError(bin);

				//printf("            bin = %i, x_cen = %E, y = %E, y_unc = %E\n", bin, x_cen, y, y_unc);

				double rel_diff = (y - norm * Gauss(x_cen, entry.si_ref)) / y_unc;

				S2 += rel_diff * rel_diff;
			}
		}

		printf("        S2 = %E\n", S2);

		int idx = g_chi_sq->GetN();
		g_chi_sq->SetPoint(idx, sh, S2);

		if (S2 < S2_best)
		{
			S2_best = S2;
			sh_best = sh;
			norm_vec_best = norm_vec;
		}
	}

	// save plots
	for (const auto &e : entries)
		e.h_th_y->Write();

	g_chi_sq->Write();

	// make final match plot
	TGraphErrors *g_hist_sh = new TGraphErrors(); g_hist_sh->SetName("g_hist_sh"); g_hist_sh->SetMarkerColor(2);
	TGraph *g_ref_fcn = new TGraph(); g_ref_fcn->SetName("g_ref_fcn"); g_ref_fcn->SetLineColor(4);
	
	for (unsigned int ei = 0; ei < entries.size(); ei++)
	{
		const Entry &entry = entries[ei];
		const EntryAux &aux = entriesAux[ei];

		for (int bin = aux.bin_min; bin <= aux.bin_max; ++bin)
		{
			int idx = g_hist_sh->GetN();
			double x = entry.h_th_y->GetBinCenter(bin) + sh_best;

			g_hist_sh->SetPoint(idx, x, entry.h_th_y->GetBinContent(bin));
			g_hist_sh->SetPointError(idx, 0., entry.h_th_y->GetBinError(bin));

			g_ref_fcn->SetPoint(idx, x, Gauss(x, entry.si_ref) * norm_vec_best[ei]);
		}
	}

	TCanvas *c = new TCanvas();
	c->SetName("cmp");
	g_hist_sh->Draw("ap");
	g_ref_fcn->Draw("l");
	c->Write();

	// restore global directory
	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// settings
	vector<string> datasets = {
		"DS1",
		"DS1_no_add_alignment",
	};
	
	// prepare output
	TFile *f_out = TFile::Open("match_th_y_scale.root", "recreate");

	// process input
	for (const auto &dataset : datasets)
	{
		TDirectory *d_dataset = f_out->mkdir(dataset.c_str());

		TFile *f_in;
		string diagonal;

		double si_ref = 124.4E-6;	// mu rad
	
		//--------------------

		diagonal = "45b_56t";	
		f_in = TFile::Open(("../" + dataset + "/distributions_" + diagonal + ".root").c_str());
		gDirectory = d_dataset->mkdir(diagonal.c_str());

		{
			vector<Entry> entries = {
				{ (TH1D *)f_in->Get("selected - angles/h_th_y_L_F"), +220E-6, +420E-6, si_ref }
			};
			DoMatch(entries, "L_F");
		}

		delete f_in;

		//--------------------

		diagonal = "45t_56b";	
		f_in = TFile::Open(("../" + dataset + "/distributions_" + diagonal + ".root").c_str());
		gDirectory = d_dataset->mkdir(diagonal.c_str());

		{
			vector<Entry> entries = {
				{ (TH1D *)f_in->Get("selected - angles/h_th_y_L_F"), -420E-6, -220E-6, si_ref }
			};
			DoMatch(entries, "L_F");
		}

		delete f_in;
	}

	// clean up
	delete f_out;

	return 0;
}
