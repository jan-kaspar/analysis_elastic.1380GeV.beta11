#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"

#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

void FitOne(const string &diagonal, TFile *f_in, const string &path, double x_from, double x_to, TF1 *ff)
{
	printf("\n\n* diagonal: %s, object: %s\n", diagonal.c_str(), path.c_str());

	ff->SetParameters(0., 0.);
	TProfile *p = (TProfile *) f_in->Get(path.c_str());
	p->Fit(ff, "", "", x_from, x_to);
	p->Write();
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// settings
	vector<string> datasets = {
		"DS1",
		"DS1_glob_al_45t_56b",
		"DS1_rel_al_only"
	};

	// prepare output
	TFile *f_out = TFile::Open("fit_th_y_diff.root", "recreate");

	// prepare processing
	TF1 *ff_pol0 = new TF1("ff_pol0", "[0]");
	TF1 *ff_pol1 = new TF1("ff_pol1", "[0] + [1]*x");

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

		FitOne(diagonal, f_in, "selected - angles/p_th_y_L_diffNF_vs_th_y_L", 220E-6, 455E-6, ff_pol1);
		FitOne(diagonal, f_in, "selected - angles/p_th_y_R_diffNF_vs_th_y_R", 220E-6, 455E-6, ff_pol1);
		FitOne(diagonal, f_in, "selected - angles/p_th_y_diffLR_vs_th_y", 220E-6, 455E-6, ff_pol0);

		delete f_in;

		//--------------------

		diagonal = "45t_56b";	
		f_in = TFile::Open(("../" + dataset + "/distributions_" + diagonal + ".root").c_str());
		gDirectory = d_dataset->mkdir(diagonal.c_str());

		FitOne(diagonal, f_in, "selected - angles/p_th_y_L_diffNF_vs_th_y_L", -455E-6, -250E-6, ff_pol1);
		FitOne(diagonal, f_in, "selected - angles/p_th_y_R_diffNF_vs_th_y_R", -455E-6, -240E-6, ff_pol1);
		FitOne(diagonal, f_in, "selected - angles/p_th_y_diffLR_vs_th_y", -350E-6, -220E-6, ff_pol0);

		delete f_in;
	}

	// clean up
	delete f_out;

	return 0;
}
