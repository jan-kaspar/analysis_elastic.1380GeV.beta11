#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"

#include <cmath>

using namespace std;

// th_x cut - for this analysis: just an additional test
double th_x_min = -50E-6;
double th_x_max = +50E-6;

//----------------------------------------------------------------------------------------------------

struct HistGroup
{
	TH1D *de_th_x, *de_th_y;
	TH1D *th_x;
	TH1D *th_y;
	TH1D *th_y_cut;
	TH2D *th_x_th_y;

	void Init()
	{
		de_th_x = new TH1D("", ";#Delta #theta_{x}", 1000, 0., 0.);
		de_th_y = new TH1D("", ";#Delta #theta_{y}", 1000, 0., 0.);

		th_x = new TH1D("", ";#theta_{x}", 500, -800E-6, 800E-6);
		th_y = new TH1D("", ";#theta_{y}", 100, 0E-6, 800E-6);
		th_y_cut = new TH1D("", ";#theta_{y}", 100, 0E-6, 500E-6);

		th_x_th_y = new TH2D("", ";#theta_{x};#theta_{y}", 60, -500E-6, 500E-6, 30, 0E-6, 500E-6);
	}

	void Fill(double thx, double thy)
	{
		if (!th_x)
			Init();
					
		th_x->Fill(thx);
		th_y->Fill(fabs(thy));

		if (thx > th_x_min && thx < th_x_max)
			th_y_cut->Fill(fabs(thy));

		th_x_th_y->Fill(thx, fabs(thy));
	}

	void FillDelta(double dthx, double dthy)
	{
		if (!de_th_x)
			Init();

		de_th_x->Fill(dthx);
		de_th_y->Fill(dthy);
	}
};

//----------------------------------------------------------------------------------------------------

TH1D* MakeSimpleRatio(TH1D *num, TH1D *den, TF1 *ff, double from, double to, bool verbose = false)
{
	TH1D *rat = new TH1D(*num);
	for (int i = 1; i <= rat->GetNbinsX(); i++)
	{
		double n = num->GetBinContent(i);
		double d = den->GetBinContent(i);

		double r = (d > 0.) ? n/d : 0.;
		double v = d * r * (1. - r);
		double e = (v >= 0. && d > 0.) ? sqrt(v) / d : 0.;

		rat->SetBinContent(i, r);
		rat->SetBinError(i, e);
	}

	rat->Fit(ff, "Q", "", from, to);
	double ep = ff->GetParameter(0);
	double ep_err = ff->GetParError(0);
	
	if (verbose)
		printf("\t\tsimple fit:  %.4f +- %.4f\n", ep, ep_err);

	return rat;
}

//----------------------------------------------------------------------------------------------------

TH1D* MakeRefinedRatio(TH1D *num, TH1D *den, TF1 *ff, double from, double to, bool verbose = false)
{
	double ep = ff->GetParameter(0);
	double ep_err = ff->GetParError(0);

	// improve error calculation: si(y)^2 = N ep (1-ep)
	TH1D *rat = new TH1D(*num);
	for (int i = 1; i <= rat->GetNbinsX(); i++)
	{
		double n = num->GetBinContent(i);
		double d = den->GetBinContent(i);

		double r = (d > 0.) ? n / d : 0.;
		double v = d * ep * (1. - ep);
		double e = (v >= 0. && d > 0.) ? sqrt(v) / d : 0.;

		rat->SetBinContent(i, r);
		rat->SetBinError(i, e);
	}
	
	rat->Fit(ff, "Q", "", from, to);
	ep = ff->GetParameter(0);
	ep_err = ff->GetParError(0);
	
	if (verbose)
		printf("\t\trefined fit: %.4f +- %.4f\n", ep, ep_err);

	return rat;
}

//----------------------------------------------------------------------------------------------------

TH2D* MakeSimpleRatio(TH2D *num, TH2D *den)
{
	TH2D *rat = new TH2D(*num);
	for (int i = 1; i <= rat->GetNbinsX(); i++)
	{
		for (int j = 1; j <= rat->GetNbinsY(); j++)
		{
			double n = num->GetBinContent(i, j);
			double d = den->GetBinContent(i, j);
	
			double r = (d > 0.) ? n/d : 0.;
			double v = d * r * (1. - r);
			double e = (v >= 0. && d > 0.) ? sqrt(v) / d : 0.;
	
			rat->SetBinContent(i, j, r);
			rat->SetBinError(i, j, e);
		}
	}

	return rat;
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal settings
	Init(argv[1]);
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	// alignment init
	for (unsigned int i = 0; i < alignmentSources.size(); ++i)
	{
		printf("\n---------- alignment source %u ----------\n", i);
		alignmentSources[i].Init();
	}
	printf("\n\n");

	// init files
	TFile *inF = new TFile((string("distill_") + argv[1] + ".root").c_str());
	TFile *outF = new TFile((string("eff3outof4_") + argv[1] + ".root").c_str(), "recreate");

	// get input data
	TTree *inT = (TTree *) inF->Get("distilled");
	EventRed ev;
	inT->SetBranchAddress("timestamp", &ev.timestamp);
	inT->SetBranchAddress("run_num", &ev.run_num);
	inT->SetBranchAddress("bunch_num", &ev.bunch_num);
	inT->SetBranchAddress("event_num", &ev.event_num);
	inT->SetBranchAddress("trigger_bits", &ev.trigger_bits);

	inT->SetBranchAddress("v_L_F", &ev.h.v_L_F); inT->SetBranchAddress("x_L_F", &ev.h.x_L_F); inT->SetBranchAddress("y_L_F", &ev.h.y_L_F);
	inT->SetBranchAddress("v_L_N", &ev.h.v_L_N); inT->SetBranchAddress("x_L_N", &ev.h.x_L_N); inT->SetBranchAddress("y_L_N", &ev.h.y_L_N);
	inT->SetBranchAddress("v_R_N", &ev.h.v_R_N); inT->SetBranchAddress("x_R_N", &ev.h.x_R_N); inT->SetBranchAddress("y_R_N", &ev.h.y_R_N);
	inT->SetBranchAddress("v_R_F", &ev.h.v_R_F); inT->SetBranchAddress("x_R_F", &ev.h.x_R_F); inT->SetBranchAddress("y_R_F", &ev.h.y_R_F);

	// tolerances (= 1 sigma of left-right difference)
	double si_de_th_x = 170E-6;
	double si_de_th_y = 30.0E-6;

	// vector of cut sigma multiples
	vector<double> n_si;
	n_si.push_back(1.);
	n_si.push_back(3.);
	n_si.push_back(5.);

	// studied RP combinations
	vector<string> rps;
	rps.push_back("L_F");
	rps.push_back("L_N");
	rps.push_back("R_N");
	rps.push_back("R_F");

	rps.push_back("L_F,R_N");
	rps.push_back("L_F,R_F");
	rps.push_back("L_N,R_N");
	rps.push_back("L_N,R_F");

	// book histograms
	vector< vector<HistGroup> > h_sel(rps.size(), vector<HistGroup>(n_si.size())), h_full(rps.size(), vector<HistGroup>(n_si.size()));

	// build histograms
	for (int ev_idx = 0; ev_idx < inT->GetEntries(); ++ev_idx)
	{
		inT->GetEntry(ev_idx);

		// skip blinded events
		if (SkipBlindedEvent(ev.event_num))
			continue;
		
		// choose the desired trigger
		if ((ev.trigger_bits & 3) == 0)	// RP trigger only
			continue;

		// remove troublesome runs
		unsigned int run = ev.run_num / 10000;
		unsigned int file = ev.run_num % 10000;
		if (SkipRun(run, file, true))
			continue;

		// select the elastic-trigger bunch(es) only
		if (SkipBunch(run, ev.bunch_num))
			continue;

		// time selection
		/*
		if (ev.timestamp < 16700 || ev.timestamp > 17300)
			continue;
			*/

		// apply fine alignment
		HitData h_al = ev.h;
		for (unsigned int i = 0; i < alignmentSources.size(); ++i)
		{
			AlignmentData alData = alignmentSources[i].Eval(ev.timestamp);
			h_al = h_al.ApplyAlignment(alData);
		}

		// loop over tested RP combinations
		for (unsigned int rpi = 0; rpi < rps.size(); rpi++)
		{
			bool L_F = (rps[rpi].find("L_F") == string::npos);
			bool L_N = (rps[rpi].find("L_N") == string::npos);
			bool R_N = (rps[rpi].find("R_N") == string::npos);
			bool R_F = (rps[rpi].find("R_F") == string::npos);

			// tracks in all selected RPs?
			if ( (L_F && !ev.h.v_L_F) || (L_N && !ev.h.v_L_N) || (R_N && !ev.h.v_R_N) || (R_F && !ev.h.v_R_F) ) 
				continue;
			
			// tracks in all 4/4 RPs?
			bool allTracks = ev.h.v_L_F && ev.h.v_L_N && ev.h.v_R_N && ev.h.v_R_F;

			// calculate elastic kinematics for the selected tracks
			double th_y_L_sel = 0., th_x_L_sel = 0., norm_L = 0.;
			if (L_F) { norm_L += 1.; th_x_L_sel += -h_al.x_L_F / env.L_x_L_F; th_y_L_sel += -h_al.y_L_F / env.L_y_L_F; }
			if (L_N) { norm_L += 1.; th_x_L_sel += -h_al.x_L_N / env.L_x_L_N; th_y_L_sel += -h_al.y_L_N / env.L_y_L_N; }
			th_x_L_sel /= norm_L; th_y_L_sel /= norm_L;
			
			double th_y_R_sel = 0., th_x_R_sel = 0., norm_R = 0.;
			if (R_F) { norm_R += 1.; th_x_R_sel += +h_al.x_R_F / env.L_x_R_F; th_y_R_sel += +h_al.y_R_F / env.L_y_R_F; }
			if (R_N) { norm_R += 1.; th_x_R_sel += +h_al.x_R_N / env.L_x_R_N; th_y_R_sel += +h_al.y_R_N / env.L_y_R_N; }
			th_x_R_sel /= norm_R; th_y_R_sel /= norm_R;

			double de_th_x_sel = th_x_R_sel - th_x_L_sel;
			double de_th_y_sel = th_y_R_sel - th_y_L_sel;

			double th_x_ref = (th_x_R_sel + th_x_L_sel) / 2.;
			double th_y_ref = (th_y_R_sel + th_y_L_sel) / 2.;
		
			// cuts on elastic kinematics at various sigma-levels
			for (unsigned int nsi = 0; nsi < n_si.size(); nsi++)
			{
				//h_sel[rpi][nsi].FillDelta(de_th_x_sel, de_th_y_sel);

				bool cut_th_x = (fabs(de_th_x_sel) < n_si[nsi] * si_de_th_x);
				bool cut_th_y = (fabs(de_th_y_sel) < n_si[nsi] * si_de_th_y);

				if (cut_th_x && cut_th_y)
				{
					h_sel[rpi][nsi].FillDelta(de_th_x_sel, de_th_y_sel);
					h_sel[rpi][nsi].Fill(th_x_ref, th_y_ref);
			
					if (allTracks)
						h_full[rpi][nsi].Fill(th_x_ref, th_y_ref);
				}
			}
		}
	}

	// process and save histograms
	TF1 *ff = new TF1("ff", "[0]");

	for (unsigned int rpi = 0; rpi < rps.size(); rpi++)
	{
		printf("\n\n------------------------- excluded RPs: %s -------------------------\n", rps[rpi].c_str());

		char buf[100];
		sprintf(buf, "excluded RPs %s", rps[rpi].c_str());
		TDirectory *rpDir = outF->mkdir(buf);

		for (unsigned int nsi = 0; nsi < n_si.size(); nsi++)
		{
			printf("* n_si: %.1f\n", n_si[nsi]);

			sprintf(buf, "n_si %.1f", n_si[nsi]);
			TDirectory *nsiDir = rpDir->mkdir(buf);

			gDirectory = nsiDir;
			h_sel[rpi][nsi].de_th_x->Write("de_th_x");
			h_sel[rpi][nsi].de_th_y->Write("de_th_y");

			TCanvas *c;

			//--------------------

			gDirectory = nsiDir->mkdir("th_x dependence");
			printf("\tth_x dependence\n");
			c = new TCanvas();
			h_sel[rpi][nsi].th_x->SetName("h_sel.th_x");
			h_sel[rpi][nsi].th_x->Sumw2();
			h_sel[rpi][nsi].th_x->SetLineColor(2);
			h_sel[rpi][nsi].th_x->Draw();
			h_full[rpi][nsi].th_x->SetName("h_full.th_x");
			h_full[rpi][nsi].th_x->Sumw2();
			h_full[rpi][nsi].th_x->SetLineColor(4);
			h_full[rpi][nsi].th_x->Draw("same");
			c->Write("th_x comparison");

			TH1D *h_simple_ratio_vs_th_x = MakeSimpleRatio(h_full[rpi][nsi].th_x, h_sel[rpi][nsi].th_x, ff, -200E-6, 200E-6, false);
			h_simple_ratio_vs_th_x->SetName("h_simple_ratio.th_x");
			h_simple_ratio_vs_th_x->Write();

			TH1D *h_refined_ratio_vs_th_x = MakeRefinedRatio(h_full[rpi][nsi].th_x, h_sel[rpi][nsi].th_x, ff, -200E-6, 200E-6, false);
			h_refined_ratio_vs_th_x->SetName("h_refined_ratio.th_x");
			h_refined_ratio_vs_th_x->Write();

			//--------------------

			gDirectory = nsiDir->mkdir("th_y dependence");
			printf("\tth_y dependence\n");
			c = new TCanvas();
			h_sel[rpi][nsi].th_y->SetName("h_sel.th_y");
			h_sel[rpi][nsi].th_y->Sumw2();
			h_sel[rpi][nsi].th_y->SetLineColor(2);
			h_sel[rpi][nsi].th_y->Draw();
			h_full[rpi][nsi].th_y->SetName("h_full.th_y");
			h_full[rpi][nsi].th_y->Sumw2();
			h_full[rpi][nsi].th_y->SetLineColor(4);
			h_full[rpi][nsi].th_y->Draw("same");
			c->Write("th_y comparison");

			TH1D *h_simple_ratio_vs_th_y = MakeSimpleRatio(h_full[rpi][nsi].th_y, h_sel[rpi][nsi].th_y, ff, anal.eff_th_y_min, 500E-6, true);
			h_simple_ratio_vs_th_y->SetName("h_simple_ratio.th_y");
			h_simple_ratio_vs_th_y->Write();

			TH1D *h_refined_ratio_vs_th_y = MakeRefinedRatio(h_full[rpi][nsi].th_y, h_sel[rpi][nsi].th_y, ff, anal.eff_th_y_min, 500E-6, true);
			h_refined_ratio_vs_th_y->SetName("h_refined_ratio.th_y");
			h_refined_ratio_vs_th_y->Write();

			//--------------------

			gDirectory = nsiDir->mkdir("th_y dependence with th_x cut");
			printf("\tth_y dependence with th_x cut\n");
			c = new TCanvas();
			h_sel[rpi][nsi].th_y_cut->SetName("h_sel.th_y");
			h_sel[rpi][nsi].th_y_cut->Sumw2();
			h_sel[rpi][nsi].th_y_cut->SetLineColor(2);
			h_sel[rpi][nsi].th_y_cut->Draw();
			h_full[rpi][nsi].th_y_cut->SetName("h_full.th_y");
			h_full[rpi][nsi].th_y_cut->Sumw2();
			h_full[rpi][nsi].th_y_cut->SetLineColor(4);
			h_full[rpi][nsi].th_y_cut->Draw("same");
			c->Write("th_y comparison");

			TH1D *h_simple_ratio_vs_th_y_cut = MakeSimpleRatio(h_full[rpi][nsi].th_y_cut, h_sel[rpi][nsi].th_y_cut, ff, anal.eff_th_y_min, 500E-6, false);
			h_simple_ratio_vs_th_y_cut->SetName("h_simple_ratio.th_y");
			h_simple_ratio_vs_th_y_cut->Write();

			TH1D *h_refined_ratio_vs_th_y_cut = MakeRefinedRatio(h_full[rpi][nsi].th_y_cut, h_sel[rpi][nsi].th_y_cut, ff, anal.eff_th_y_min, 500E-6, false);
			h_refined_ratio_vs_th_y_cut->SetName("h_refined_ratio.th_y");
			h_refined_ratio_vs_th_y_cut->Write();

			//--------------------

			gDirectory = nsiDir->mkdir("th_x, th_y dependence");
			h_sel[rpi][nsi].th_x_th_y->SetName("h_sel.th_x_th_y");
			h_sel[rpi][nsi].th_x_th_y->Write();
			h_full[rpi][nsi].th_x_th_y->SetName("h_full.th_x_th_y");
			h_full[rpi][nsi].th_x_th_y->Write();

			TH2D *h_simple_ratio_vs_th_x_th_y = MakeSimpleRatio(h_full[rpi][nsi].th_x_th_y, h_sel[rpi][nsi].th_x_th_y);
			h_simple_ratio_vs_th_x_th_y->SetName("h_simple_ratio_vs_th_x_th_y");
			h_simple_ratio_vs_th_x_th_y->Write();
		}
	}
	
	return 0;
}
