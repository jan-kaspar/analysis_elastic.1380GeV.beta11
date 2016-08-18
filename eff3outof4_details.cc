#include "input_files.h"

#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

#include <vector>
#include <map>
#include <string>
#include <cmath>

#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"
#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/RawDataFormats.h"
#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/RPRootTrackInfo.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

struct RPStruct
{
	RPRootDumpDigiInfo *digi;
	RPRootDumpPatternInfo *pat;
	RPRootDumpTrackInfo *tr;
	vector<RPRootDumpTrackInfo> *mtr;

	RPStruct() : digi(NULL), pat(NULL), tr(NULL) {}

	void AssignBranches(TChain *ch, unsigned int id) {
		char buf[100];

		digi = new RPRootDumpDigiInfo;
		sprintf(buf, "digi_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "digi_rp_%u.", id); ch->SetBranchAddress(buf, &digi);

		pat = new RPRootDumpPatternInfo();
		sprintf(buf, "nonpar_patterns_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "nonpar_patterns_rp_%u.", id); ch->SetBranchAddress(buf, &pat);

		tr = new RPRootDumpTrackInfo();
		sprintf(buf, "track_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch->SetBranchAddress(buf, &tr);

		//mtr = new vector<RPRootDumpTrackInfo>();
		//sprintf(buf, "multi_track_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		//sprintf(buf, "multi_track_rp_%u", id); ch->SetBranchAddress(buf, &mtr);
	}

	/*
	double x, y;		// track hit including alignment
	double th_x, th_y;

	void SimpleReco(double a, double b, double c, double Lx, double Ly) {
		x = tr->x - a * tr->y - b;
	   	y = tr->y - c;
		th_x = x / Lx;
		th_y = y / Ly;
	}
	*/
};

//----------------------------------------------------------------------------------------------------

struct DiagStruct
{
	RPStruct L_F, L_N, R_N, R_F;

	void AssignBranches(TChain *ch, unsigned int lf, unsigned int ln, unsigned int rn, unsigned int rf) {
		L_F.AssignBranches(ch, lf);
		L_N.AssignBranches(ch, ln);
		R_N.AssignBranches(ch, rn);
		R_F.AssignBranches(ch, rf);
	}
};

//----------------------------------------------------------------------------------------------------

struct HistGroup
{
	TH1D *h_y, *h_th_x, *h_th_y;
	TH2D *h_th_x_th_y;

	HistGroup() : h_y(NULL), h_th_y(NULL) {}

	void Init()
	{
		h_y = new TH1D("", ";y   (mm)", 100, -40., +40.);

		h_th_x = new TH1D("", ";#theta_{x}   (#murad)", 300, -600., +600.);
		h_th_y = new TH1D("", ";#theta_{y}   (#murad)", 300, -600., +600.);

		h_th_x_th_y = new TH2D("", ";#theta_{x}   (#murad);#theta_{y}   (#murad)", 30, -600, 600, 30, -600, 600);
	}

	void Fill(double y, double th_x, double th_y)
	{
		if (!h_y)
			Init();

		h_y->Fill(y);

		h_th_x->Fill(th_x * 1E6);
		h_th_y->Fill(th_y * 1E6);
		h_th_x_th_y->Fill(th_x * 1E6, th_y * 1E6);
	}

	void Write()
	{
		h_y->Write();

		h_th_x->Write();
		h_th_y->Write();
		h_th_x_th_y->Write();
	}
};

//----------------------------------------------------------------------------------------------------

// map: element set, condition -> group of histograms
typedef map<string, map<string, HistGroup> > CounterMap;

//----------------------------------------------------------------------------------------------------

/**
NON-PARALLEL pattern recognition run with these parameters:
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 3
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.threshold = 2.99
**/

bool RPTooFullU(const RPStruct &rp)
{
	// count planes that could have (in principle) contributed to pattern search
	unsigned N = 0;
	for (unsigned int i = 1; i < 10; i += 2)
		if (rp.digi->numberOfClusters[i] <= 5)
			N++;

	return (N < 3);
}

//----------------------------------------------------------------------------------------------------

bool RPTooFullV(const RPStruct &rp)
{
	unsigned N = 0;
	for (unsigned int i = 0; i < 10; i += 2)
		if (rp.digi->numberOfClusters[i] <= 5)
			N++;

	return (N < 3);
}

//----------------------------------------------------------------------------------------------------

void AnalyzeOnePot(const string &pot, const DiagStruct &dgn, const vector<AlignmentData> &alData, CounterMap &c)
{
	// which pots are selected for track/event definition
	bool sel_L_F = (pot.find("L_F") == string::npos);
	bool sel_L_N = (pot.find("L_N") == string::npos);
	bool sel_R_N = (pot.find("R_N") == string::npos);
	bool sel_R_F = (pot.find("R_F") == string::npos);

	// do all selected pots have tracks?
	bool skip = (sel_L_F && !dgn.L_F.tr->valid) || (sel_L_N && !dgn.L_N.tr->valid)
		|| (sel_R_N && !dgn.R_N.tr->valid) || (sel_R_F && !dgn.R_F.tr->valid);

	if (skip)
		return;

	// apply alignment
	HitData h_raw;
	h_raw.x_L_F = dgn.L_F.tr->x; h_raw.y_L_F = dgn.L_F.tr->y;
	h_raw.x_L_N = dgn.L_N.tr->x; h_raw.y_L_N = dgn.L_N.tr->y;
	h_raw.x_R_N = dgn.R_N.tr->x; h_raw.y_R_N = dgn.R_N.tr->y;
	h_raw.x_R_F = dgn.R_F.tr->x; h_raw.y_R_F = dgn.R_F.tr->y;
	
	HitData h_al = h_raw;
	for (unsigned int i = 0; i < alData.size(); i++)
		h_al = h_al.ApplyAlignment(alData[i]);

	// calculate elastic kinematics for the selected tracks
	double th_y_L_sel = 0., th_x_L_sel = 0., norm_L = 0.;
	if (sel_L_F) { norm_L += 1.; th_x_L_sel += -h_al.x_L_F / env.L_x_L_F; th_y_L_sel += -h_al.y_L_F / env.L_y_L_F; }
	if (sel_L_N) { norm_L += 1.; th_x_L_sel += -h_al.x_L_N / env.L_x_L_N; th_y_L_sel += -h_al.y_L_N / env.L_y_L_N; }
	th_x_L_sel /= norm_L; th_y_L_sel /= norm_L;
	
	double th_y_R_sel = 0., th_x_R_sel = 0., norm_R = 0.;
	if (sel_R_F) { norm_R += 1.; th_x_R_sel += +h_al.x_R_F / env.L_x_R_F; th_y_R_sel += +h_al.y_R_F / env.L_y_R_F; }
	if (sel_R_N) { norm_R += 1.; th_x_R_sel += +h_al.x_R_N / env.L_x_R_N; th_y_R_sel += +h_al.y_R_N / env.L_y_R_N; }
	th_x_R_sel /= norm_R; th_y_R_sel /= norm_R;

	double th_x_sel = (th_x_R_sel + th_x_L_sel) / 2.;
	double th_y_sel = (th_y_R_sel + th_y_L_sel) / 2.;

	double de_th_x_sel = th_x_R_sel - th_x_L_sel;
	double de_th_y_sel = th_y_R_sel - th_y_L_sel;

	// can this be elastic event ?
	double n_si = 3.;
	double si_de_th_x = 170E-6;
	double si_de_th_y = 30E-6;

	bool cut_th_x = (fabs(de_th_x_sel) < n_si * si_de_th_x);
	bool cut_th_y = (fabs(de_th_y_sel) < n_si * si_de_th_y);

	if (!cut_th_x || !cut_th_y)
		return;

	// pot under test
	RPStruct rp_test;
	double y_test = 0., th_x_test = 0., th_y_test = 0.;
	if (pot.compare("L_F") == 0) { rp_test = dgn.L_F; y_test = - h_al.y_L_F; th_x_test = - h_al.x_L_F / env.L_x_L_F; th_y_test = - h_al.y_L_F / env.L_y_L_F; }
	if (pot.compare("L_N") == 0) { rp_test = dgn.L_N; y_test = - h_al.y_L_N; th_x_test = - h_al.x_L_N / env.L_x_L_N; th_y_test = - h_al.y_L_N / env.L_y_L_N; }
	if (pot.compare("R_N") == 0) { rp_test = dgn.R_N; y_test = + h_al.y_R_N; th_x_test = + h_al.x_R_N / env.L_x_R_N; th_y_test = + h_al.y_R_N / env.L_y_R_N; }
	if (pot.compare("R_F") == 0) { rp_test = dgn.R_F; y_test = + h_al.y_R_F; th_x_test = + h_al.x_R_F / env.L_x_R_F; th_y_test = + h_al.y_R_F / env.L_y_R_F; }

	// reference quantities
	double th_x_ref = th_x_sel;
	double th_y_ref = th_y_sel;
	double y_ref = y_test;

	// what happens with the test pot
	bool rp_test_pl_insuff = (rp_test.digi->uPlanesOn < 3 || rp_test.digi->vPlanesOn < 3);
	bool rp_test_pl_suff = (rp_test.digi->uPlanesOn >= 3 && rp_test.digi->vPlanesOn >= 3);
	bool rp_test_pl_suff_no_track = rp_test_pl_suff && !rp_test.tr->valid;
	bool rp_test_pl_too_full_u = RPTooFullU(rp_test), rp_test_pl_too_full_v = RPTooFullV(rp_test);
	//bool rp_test_pat_suff = (rp_test.pat->u.size() > 0 || rp_test.pat->v.size() > 0) || (rp_test_pl_too_full_u || rp_test_pl_too_full_v);
	bool rp_test_pat_more = (rp_test.pat->u.size() > 1 && rp_test.pat->v.size() > 1) || (rp_test_pl_too_full_u || rp_test_pl_too_full_v);
	bool rp_test_track = rp_test.tr->valid;
	bool rp_test_track_compatible = rp_test_track
		//&& (fabs(th_x_test - th_x_ref) < n_si * si_de_th_x)
		&& (fabs(th_y_test - th_y_ref) < n_si * si_de_th_y);

	c[pot]["anything"].Fill(y_ref, th_x_ref, th_y_ref);
	
	if (rp_test_pl_insuff)
		c[pot]["pl_insuff"].Fill(y_ref, th_x_ref, th_y_ref);
	if (rp_test_pl_suff_no_track)
		c[pot]["pl_suff_no_track"].Fill(y_ref, th_x_ref, th_y_ref);
	if (rp_test_pat_more)
		c[pot]["pat_more"].Fill(y_ref, th_x_ref, th_y_ref);
	if (rp_test_track)
		c[pot]["track"].Fill(y_ref, th_x_ref, th_y_ref);
	if (rp_test_track_compatible)
		c[pot]["track_compatible"].Fill(y_ref, th_x_ref, th_y_ref);
}

//----------------------------------------------------------------------------------------------------

void AnalyzeDiagonal(const DiagStruct &dgn, const vector<AlignmentData> &al, CounterMap &c)
{
	AnalyzeOnePot("L_F", dgn, al, c);
	AnalyzeOnePot("L_N", dgn, al, c);
	AnalyzeOnePot("R_N", dgn, al, c);
	AnalyzeOnePot("R_F", dgn, al, c);
}

//----------------------------------------------------------------------------------------------------

TH1D* MakeSimpleRatio(TH1D *num, TH1D *den)
{
	TH1D *rat = new TH1D(*num);

	for (int i = 1; i <= rat->GetNbinsX(); i++) {
		double n = num->GetBinContent(i);
		double d = den->GetBinContent(i);

		double r = (d > 0.) ? n/d : 0.;
		double v = d * r * (1. - r);
		double e = (v >= 0. && d > 0.) ? sqrt(v) / d : 0.;

		rat->SetBinContent(i, r);
		rat->SetBinError(i, e);
	}

	return rat;
}

//----------------------------------------------------------------------------------------------------

TH2D* MakeSimpleRatio(TH2D *num, TH2D *den)
{
	TH2D *rat = new TH2D(*num);
	for (int i = 1; i <= rat->GetNbinsX(); i++) {
		for (int j = 1; j <= rat->GetNbinsY(); j++) {
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
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal
	Init(argv[1]);
	if (diagonal != dCombined)
		return rcIncompatibleDiagonal;
	
	// alignment init
	for (unsigned int i = 0; i < alignmentSources.size(); ++i)
	{
		printf("\n---------- alignment source %u ----------\n", i);
		alignmentSources[i].Init();
	}
	printf("\n\n");

	// get input
	InitInputFiles();
	TChain *ch = new TChain("TotemNtuple");
	for (unsigned int i = 0; i < input_files.size(); i++) {
		ch->Add(input_files[i].c_str());
		printf("+ %s\n", input_files[i].c_str());
	}
	
	// prepare output
	TFile *outF = new TFile("eff3outof4_details.root", "recreate");

	// select and link input branches
	ch->SetBranchStatus("*", 0);

	EventMetaData *metaData = new EventMetaData();
	ch->SetBranchStatus("event_info.*", 1);
	ch->SetBranchAddress("event_info.", &metaData);

	TriggerData *triggerData = new TriggerData();
	ch->SetBranchStatus("trigger_data.*", 1);
	ch->SetBranchAddress("trigger_data.", &triggerData);

	DiagStruct diag_45b, diag_45t;
	diag_45b.AssignBranches(ch, 25, 21, 120, 124);
	diag_45t.AssignBranches(ch, 24, 20, 121, 125);

	// prepare counters and histograms
	map<string, CounterMap> counters;	// map: diagonal label -> CounterMap

	signed int offset = 0;
	unsigned int prev_run = 0;

	// loop over events
	for (unsigned int ev = 0; ev < ch->GetEntries(); ev++)
	{
		// load all data for event ev
		ch->GetEvent(ev);
		
		signed int bunch_ei = metaData->optoRx_BX[0] - 264;
		if (bunch_ei < 0)
			bunch_ei += 3564;

		// reset offset for new run
		unsigned int run = metaData->run_no / 10000;
		unsigned int file = metaData->run_no % 10000;
		unsigned int event = metaData->event_no;
		if (run != prev_run)
		{
			prev_run = run;
			offset = 0;
			printf(">> new run (%u) found, offset reset\n", run);
		}
		
		// load all data for event ev-offset
		if (offset != 0)
			ch->GetEvent(ev-offset);

		TriggerData td = *triggerData;
		
		// mismatch check and correction
		signed int bunch_tr = td.bunch_num;

		if (bunch_ei != bunch_tr)
		{
			printf(">> mismatch found (event_info: run %u, file %u, event %u): bunch_ei = %i, bunch_tr_corr = %i\n",
				run, file, event, bunch_ei, bunch_tr);

			printf("\toffset %u", offset);
			offset++;
			printf(" => %u\n", offset);

			printf("\tskipping event\n");
			continue;
		}

		// reload all data for event ev
		if (offset != 0)
			ch->GetEvent(ev);

		// select RP trigger only
		if ((td.input_status_bits & 3) == 0)
			continue;

		// skip blinded events
		if (SkipBlindedEvent(td.event_num))
			continue;

		// skip troublesome runs
		if (SkipRun(run, file, false))
			continue;

		// select bunches
		if (SkipBunch(run, td.bunch_num))
			continue;

		// get fine alignment
		vector<AlignmentData> alData;
		for (unsigned int i = 0; i < alignmentSources.size(); ++i)
			alData.push_back(alignmentSources[i].Eval(metaData->timestamp - timestamp0));

		// run analysis
		AnalyzeDiagonal(diag_45b, alData, counters["45b_56t"]);
		AnalyzeDiagonal(diag_45t, alData, counters["45t_56b"]);
	}

	// save results
	for (map<string, CounterMap>::iterator dgni = counters.begin(); dgni != counters.end(); ++dgni) {
		TDirectory *dgnDir = outF->mkdir(dgni->first.c_str());
		
		for (CounterMap::iterator oi = dgni->second.begin(); oi != dgni->second.end(); ++oi) {
			TDirectory *objDir = dgnDir->mkdir(oi->first.c_str());

			TH1D *h_th_x_tot = new TH1D(*oi->second["anything"].h_th_x);
			TH1D *h_th_y_tot = new TH1D(*oi->second["anything"].h_th_y);
			TH2D *h_th_x_th_y_tot = new TH2D(*oi->second["anything"].h_th_x_th_y);

			for (map<string, HistGroup>::iterator ci = oi->second.begin(); ci != oi->second.end(); ++ci) {
				TDirectory *condDir = objDir->mkdir(ci->first.c_str());
				gDirectory = condDir;

				char buf[100];

				//--------------------

				sprintf(buf, "th_x");
				ci->second.h_th_x->SetName(buf);
				ci->second.h_th_x->Write();
				
				TH1D *h_rat_th_x = MakeSimpleRatio(ci->second.h_th_x, h_th_x_tot);
				sprintf(buf, "th_x : rel");
				h_rat_th_x->SetName(buf);
				h_rat_th_x->Write();

				//--------------------

				sprintf(buf, "th_y");
				ci->second.h_th_y->SetName(buf);
				ci->second.h_th_y->Write();
				
				TH1D *h_rat_th_y = MakeSimpleRatio(ci->second.h_th_y, h_th_y_tot);
				sprintf(buf, "th_y : rel");
				h_rat_th_y->SetName(buf);
				h_rat_th_y->Write();

				//--------------------

				sprintf(buf, "th_x, th_y");
				ci->second.h_th_x_th_y->SetName(buf);
				ci->second.h_th_x_th_y->Write();
				
				TH2D *h_rat_th_x_th_y = MakeSimpleRatio(ci->second.h_th_x_th_y, h_th_x_th_y_tot);
				sprintf(buf, "th_x, th_y : rel");
				h_rat_th_x_th_y->SetName(buf);
				h_rat_th_x_th_y->Write();
			}

			delete h_th_x_tot;
			delete h_th_y_tot;
			delete h_th_x_th_y_tot;
		}
	}

	delete outF;
	return 0;
}
