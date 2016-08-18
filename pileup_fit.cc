#include "../common_definitions.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include <map>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

void DoFit(TFile *inF, const string &dir, const string &label)
{
	TGraphErrors *g_val = (TGraphErrors *) inF->Get((dir + "/rel").c_str());
	TGraphErrors *g_run = (TGraphErrors *) inF->Get((dir + "/runs").c_str());

	g_run->SetName((label+".runs").c_str());
	g_run->Write();
	
	g_val->SetName((label+".src").c_str());
	g_val->Write();

	// build graph per run
	map<unsigned int, TGraphErrors *> gMap;
	for (int pi = 0; pi < g_run->GetN(); pi++) {
		// load data from input graphs
		double t, v, runD;
		g_run->GetPoint(pi, t, runD);
		g_val->GetPoint(pi, t, v);
		double ve = g_val->GetErrorY(pi);
		double te = g_val->GetErrorX(pi);
		unsigned int run = (unsigned int) runD;

		// crate new graph, if needed
		map<unsigned int, TGraphErrors *>::iterator it = gMap.find(run);
		if (it == gMap.end())
			it = gMap.insert(pair<unsigned int, TGraphErrors *>(run, new TGraphErrors)).first;

		// fill in the input data
		int idx = it->second->GetN();
		it->second->SetPoint(idx, t, v);
		it->second->SetPointError(idx, te, ve);
	}

	// fit for each run and compile global graph
	TF1 *pol0 = new TF1("pol0", "[0]");
	TF1 *pol1 = new TF1("pol1", "[0] + [1]*x");

	TGraph *g_gl = new TGraphErrors();
	g_gl->SetName(label.c_str());
	g_gl->SetLineColor(2);
	for (map<unsigned int, TGraphErrors *>::iterator pit = gMap.begin(); pit != gMap.end(); ++pit) {
		//printf("* %s / %u\n", label.c_str(), pit->first);
		TF1 *ff = (pit->second->GetN() > 1) ? pol1 : pol0;

		// fit the graph from one run
		pit->second->Fit(ff, "WQ");
		pit->second->Fit(ff, "Q");

		// determine boundaries of the run
		double t_min=0., t_max=0.;
		for (int i = 0; i < pit->second->GetN(); i++) {
			double x, y;
			pit->second->GetPoint(i, x, y);
			double ex = pit->second->GetErrorX(i);

			t_min = (t_min == 0.) ? x-ex : min(t_min, x-ex);
			t_max = (t_max == 0.) ? x+ex : max(t_max, x+ex);
		}

		// add the fit points at the boundaries to the global graph
		g_gl->SetPoint(g_gl->GetN(), t_min, ff->Eval(t_min));
		g_gl->SetPoint(g_gl->GetN(), t_max, ff->Eval(t_max));

		delete pit->second;
	}
	g_gl->Write();

	TCanvas *c = new TCanvas((label+".cmp").c_str());
	g_val->Draw("ap");
	g_gl->Draw("l");
	c->Write();
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

	// get input
	TFile *inF = new TFile("pileup_combined.root");
	
	// prepare output
	TFile *outF = new TFile((string("pileup_fit_") + argv[1] + ".root").c_str(), "recreate");

	gDirectory = outF->mkdir("45b_56t");
	DoFit(inF, "45b/dgn/pat_suff && pat_suff, L || R", "dgn");

	gDirectory = outF->mkdir("45t_56b");
	DoFit(inF, "45t/dgn/pat_suff && pat_suff, L || R", "dgn");

	delete outF;
	return 0;
}
