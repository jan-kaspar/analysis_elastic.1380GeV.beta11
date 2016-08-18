#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TH1D.h"
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

// map: element set, condition -> count 
typedef map<string, map<string, unsigned long > > CounterMap;

//----------------------------------------------------------------------------------------------------

#if 0
/**
NON-PARALLEL pattern recognition run with these parameters:
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 3
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.threshold = 2.99
**/

#endif


//----------------------------------------------------------------------------------------------------

struct RPStruct
{
	string name;

	RPRootDumpDigiInfo *digi;
	RPRootDumpPatternInfo *pat;
	RPRootDumpTrackInfo *tr;
	
	RPRootDumpTrackInfo *tr_id;

	bool prot;
	bool pl_suff;
	bool pl_too_full_u, pl_too_full_v;
	bool pat_suff;
	bool pat_more;
	//bool tr_any;
	bool tr_val;

	void AssignBranches(const string &_n, TChain *ch_full, TChain *ch_ideal, unsigned int id)
	{
		name = _n;

		char buf[100];

		digi = new RPRootDumpDigiInfo;
		sprintf(buf, "digi_rp_%u.*", id); ch_full->SetBranchStatus(buf, 1);
		sprintf(buf, "digi_rp_%u.", id); ch_full->SetBranchAddress(buf, &digi);

		pat = new RPRootDumpPatternInfo();
		sprintf(buf, "nonpar_patterns_rp_%u.*", id); ch_full->SetBranchStatus(buf, 1);
		sprintf(buf, "nonpar_patterns_rp_%u.", id); ch_full->SetBranchAddress(buf, &pat);

		tr = new RPRootDumpTrackInfo();
		sprintf(buf, "track_rp_%u.*", id); ch_full->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch_full->SetBranchAddress(buf, &tr);

		tr_id = new RPRootDumpTrackInfo();
		sprintf(buf, "track_rp_%u.*", id); ch_ideal->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch_ideal->SetBranchAddress(buf, &tr_id);
	}
	
	bool RPTooFullU()
	{
		// count planes that could have (in principle) contributed to pattern search
		unsigned N = 0;
		for (unsigned int i = 1; i < 10; i += 2)
			if (digi->numberOfClusters[i] <= 5)
				N++;
	
		return (N < 3);
	}
	
	bool RPTooFullV()
	{
		unsigned N = 0;
		for (unsigned int i = 0; i < 10; i += 2)
			if (digi->numberOfClusters[i] <= 5)
				N++;
	
		return (N < 3);
	}
	void Analyze(CounterMap &c);
};

//----------------------------------------------------------------------------------------------------

void RPStruct::Analyze(CounterMap &c)
{
	prot = tr_id->valid;

	pl_suff = (digi->uPlanesOn >= 3 && digi->vPlanesOn >= 3);
	pl_too_full_u = RPTooFullU();
	pl_too_full_v = RPTooFullV();
	pat_suff = (pat->u.size() > 0 || pat->v.size() > 0) || (pl_too_full_u || pl_too_full_v);
	pat_more = (pat->u.size() > 1 && pat->v.size() > 1) || (pl_too_full_u || pl_too_full_v);
	//tr_any = (tr->valid || mtr->size() > 0);

	tr_val = tr->valid;

	if (prot) c[name]["prot"]++;

	if (prot & pat_suff) c[name]["prot & pat_suff"]++;
	if (prot & !pat_suff) c[name]["prot & !pat_suff"]++;
	if (!prot & pat_suff) c[name]["!prot & pat_suff"]++;
	if (!prot & !pat_suff) c[name]["!prot & !pat_suff"]++;

	if (prot & tr_val) c[name]["prot & tr_val"]++;
	if (prot & !tr_val) c[name]["prot & !tr_val"]++;
	if (!prot & tr_val) c[name]["!prot & tr_val"]++;
	if (!prot & !tr_val) c[name]["!prot & !tr_val"]++;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct StationStruct
{
	string name;

	RPStruct N_T, N_B, F_T, F_B;

	map<string, TH1D *> hists_th_y;

	void AssignBranches(const string &_n, TChain *ch_full, TChain *ch_ideal,
			unsigned int nt, unsigned int nb, unsigned int ft, unsigned int fb)
	{
		name = _n;

		N_T.AssignBranches("N_T", ch_full, ch_ideal, nt);
		N_B.AssignBranches("N_B", ch_full, ch_ideal, nb);
		F_T.AssignBranches("F_T", ch_full, ch_ideal, ft);
		F_B.AssignBranches("F_B", ch_full, ch_ideal, fb);
	}

	void FillThYHist(const string &desc, double th_y)
	{
		map<string, TH1D *>::iterator it = hists_th_y.find(desc);
		if (it == hists_th_y.end())
			it = hists_th_y.insert({desc, new TH1D("", ";th_y", 250, -500E-6, +500E-6)}).first;

		it->second->Fill(th_y);
	}

	void MakeThYHistRatio(const string &desc_num, const string &desc_den)
	{
		TH1D *h = new TH1D(* hists_th_y[desc_num]);
		h->Sumw2();

		TH1D *h_d = hists_th_y[desc_den];

		for (int bi = 1; bi <= h->GetNbinsX(); bi++)
		{
			double n = h->GetBinContent(bi);
			//double n_u = h->GetBinError(bi);
			double d = h_d->GetBinContent(bi);

			double r = (d > 0) ? n / d : 0.;
			double r_u = (d > 0) ? sqrt(n * (1. - n/d)) / d : 0.;

			h->SetBinContent(bi, r);
			h->SetBinError(bi, r_u);
		}

		hists_th_y[desc_num + " / " + desc_den] = h;
	}

	void Analyze(unsigned long ev, double th_y, CounterMap &c);

	void MakeRatios();

	void Write()
	{
		TDirectory *topDir = gDirectory;
		gDirectory = topDir->mkdir(name.c_str());

		for (map<string, TH1D *>::iterator it = hists_th_y.begin(); it != hists_th_y.end(); ++it)
		{
			it->second->Write(it->first.c_str());
		}

		gDirectory = topDir;
	}
};

//----------------------------------------------------------------------------------------------------

void StationStruct::Analyze(unsigned long ev, double th_y, CounterMap &c)
{
	c["total"]["total"]++;

	N_T.Analyze(c);
	N_B.Analyze(c);
	F_T.Analyze(c);
	F_B.Analyze(c);

	if (N_T.prot && F_T.prot) c["N_T,F_T"]["prot & prot"]++;
	if (N_T.prot && F_T.prot && N_T.tr_val && F_T.tr_val) c["N_T,F_T"]["prot & prot & tr_val & tr_val"]++;
	if (N_T.prot && F_T.prot && N_T.tr_val && !F_T.tr_val) c["N_T,F_T"]["prot & prot & tr_val & !tr_val"]++;
	if (N_T.prot && F_T.prot && !N_T.tr_val && F_T.tr_val) c["N_T,F_T"]["prot & prot & !tr_val & tr_val"]++;
	if (N_T.prot && F_T.prot && !N_T.tr_val && !F_T.tr_val) c["N_T,F_T"]["prot & prot & !tr_val & !tr_val"]++;

	if (N_B.prot && F_B.prot) c["N_B,F_B"]["prot & prot"]++;
	if (N_B.prot && F_B.prot && N_B.tr_val && F_B.tr_val) c["N_B,F_B"]["prot & prot & tr_val & tr_val"]++;
	if (N_B.prot && F_B.prot && N_B.tr_val && !F_B.tr_val) c["N_B,F_B"]["prot & prot & tr_val & !tr_val"]++;
	if (N_B.prot && F_B.prot && !N_B.tr_val && F_B.tr_val) c["N_B,F_B"]["prot & prot & !tr_val & tr_val"]++;
	if (N_B.prot && F_B.prot && !N_B.tr_val && !F_B.tr_val) c["N_B,F_B"]["prot & prot & !tr_val & !tr_val"]++;
	
	if (N_T.prot && F_T.prot && N_B.tr_val) c["N_T,F_T,N_B"]["prot & prot & opp-tr_val"]++;
	if (N_T.prot && F_T.prot && F_B.tr_val) c["N_T,F_T,F_B"]["prot & prot & opp-tr_val"]++;
	if (N_T.prot && F_T.prot && N_B.tr_val && F_B.tr_val) c["N_T,F_T,N_B,F_B"]["prot & prot & opp-tr_val & opp-tr_val"]++;
	
	if (N_T.prot && F_T.prot && N_B.pl_suff) c["N_T,F_T,N_B"]["prot & prot & opp-pl_suff"]++;
	if (N_T.prot && F_T.prot && F_B.pl_suff) c["N_T,F_T,F_B"]["prot & prot & opp-pl_suff"]++;
	if (N_T.prot && F_T.prot && N_B.pl_suff && F_B.pl_suff) c["N_T,F_T,N_B,F_B"]["prot & prot & opp-pl_suff & opp-pl_suff"]++;
	
	if (N_B.prot && F_B.prot && N_T.tr_val) c["N_B,F_B,N_T"]["prot & prot & opp-tr_val"]++;
	if (N_B.prot && F_B.prot && F_T.tr_val) c["N_B,F_B,F_T"]["prot & prot & opp-tr_val"]++;
	if (N_B.prot && F_B.prot && N_T.tr_val && F_T.tr_val) c["N_B,F_B,N_T,F_T"]["prot & prot & opp-tr_val & opp-tr_val"]++;
	
	if (N_B.prot && F_B.prot && N_T.pl_suff) c["N_B,F_B,N_T"]["prot & prot & opp-pl_suff"]++;
	if (N_B.prot && F_B.prot && F_T.pl_suff) c["N_B,F_B,F_T"]["prot & prot & opp-pl_suff"]++;
	if (N_B.prot && F_B.prot && N_T.pl_suff && F_T.pl_suff) c["N_B,F_B,N_T,F_T"]["prot & prot & opp-pl_suff & opp-pl_suff"]++;
	
	if (N_T.prot && F_T.prot) FillThYHist("N_T.prot & N_F.prot", th_y);
	if (N_T.prot && F_T.prot && N_T.tr_val && F_T.tr_val) FillThYHist("N_T.prot & F_T.prot & N_T.tr_val & F_T.tr_val", th_y);
	if (N_T.prot && F_T.prot && N_T.tr_val && !F_T.tr_val) FillThYHist("N_T.prot & F_T.prot & N_T.tr_val & !F_T.tr_val", th_y);
	if (N_T.prot && F_T.prot && !N_T.tr_val && F_T.tr_val) FillThYHist("N_T.prot & F_T.prot & !N_T.tr_val & F_T.tr_val", th_y);
	if (N_T.prot && F_T.prot && !N_T.tr_val && !F_T.tr_val) FillThYHist("N_T.prot & F_T.prot & !N_T.tr_val & !F_T.tr_val", th_y);

	if (N_B.prot && F_B.prot) FillThYHist("N_B.prot & N_F.prot", th_y);
	if (N_B.prot && F_B.prot && N_B.tr_val && F_B.tr_val) FillThYHist("N_B.prot & F_B.prot & N_B.tr_val & F_B.tr_val", th_y);
	if (N_B.prot && F_B.prot && N_B.tr_val && !F_B.tr_val) FillThYHist("N_B.prot & F_B.prot & N_B.tr_val & !F_B.tr_val", th_y);
	if (N_B.prot && F_B.prot && !N_B.tr_val && F_B.tr_val) FillThYHist("N_B.prot & F_B.prot & !N_B.tr_val & F_B.tr_val", th_y);
	if (N_B.prot && F_B.prot && !N_B.tr_val && !F_B.tr_val) FillThYHist("N_B.prot & F_B.prot & !N_B.tr_val & !F_B.tr_val", th_y);

	// event print
	/*
	if (!N_T.prot && N_T.tr_val)
		printf("\tevent_selection.insert(%lu);\t// %s, N_T\n", ev, name.c_str());

	if (!F_T.prot && F_T.tr_val)
		printf("\tevent_selection.insert(%lu);\t// %s, F_T\n", ev, name.c_str());
	*/
	
	/*
	if (N_T.prot && F_T.prot && N_T.tr_val && !F_T.tr_val)
		printf("\tevent_selection.insert({%lu, \"%s, N_T & F_T\"});\n", ev, name.c_str());
	if (N_B.prot && F_B.prot && N_B.tr_val && !F_B.tr_val)
		printf("\tevent_selection.insert({%lu, \"%s, N_B & F_B\"});\n", ev, name.c_str());
	*/

	/*
	if (N_T.prot && F_T.prot && !N_T.tr_val && F_T.tr_val)
		printf("\tevent_selection.insert(%lu);\t// %s, N_T & F_T\n", ev, name.c_str());
	if (N_B.prot && F_B.prot && !N_B.tr_val && F_B.tr_val)
		printf("\tevent_selection.insert(%lu);\t// %s, N_B & F_B\n", ev, name.c_str());
	*/

	/*
	if (N_T.prot && F_T.prot && !N_T.tr_val && !F_T.tr_val)
		printf("\tevent_selection.insert({%lu, \"%s, N_T & F_T\"});\n", ev, name.c_str());
	if (N_B.prot && F_B.prot && !N_B.tr_val && !F_B.tr_val)
		printf("\tevent_selection.insert({%lu, \"%s, N_B & F_B\"});\n", ev, name.c_str());
	*/
}

//----------------------------------------------------------------------------------------------------

void StationStruct::MakeRatios()
{
	MakeThYHistRatio("N_T.prot & F_T.prot & N_T.tr_val & F_T.tr_val", "N_T.prot & N_F.prot");
	MakeThYHistRatio("N_T.prot & F_T.prot & N_T.tr_val & !F_T.tr_val", "N_T.prot & N_F.prot");
	MakeThYHistRatio("N_T.prot & F_T.prot & !N_T.tr_val & F_T.tr_val", "N_T.prot & N_F.prot");
	MakeThYHistRatio("N_T.prot & F_T.prot & !N_T.tr_val & !F_T.tr_val", "N_T.prot & N_F.prot");

	MakeThYHistRatio("N_B.prot & F_B.prot & N_B.tr_val & F_B.tr_val", "N_B.prot & N_F.prot");
	MakeThYHistRatio("N_B.prot & F_B.prot & N_B.tr_val & !F_B.tr_val", "N_B.prot & N_F.prot");
	MakeThYHistRatio("N_B.prot & F_B.prot & !N_B.tr_val & F_B.tr_val", "N_B.prot & N_F.prot");
	MakeThYHistRatio("N_B.prot & F_B.prot & !N_B.tr_val & !F_B.tr_val", "N_B.prot & N_F.prot");
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	TChain *ch_full = new TChain("TotemNtuple");
	//ch_full->Add("ntuple_full_1E5.root");
	ch_full->Add("ntuple_full_2E5_seed1.root");
	ch_full->Add("ntuple_full_2E5_seed2.root");
	ch_full->Add("ntuple_full_2E5_seed3.root");

	TChain *ch_ideal = new TChain("TotemNtuple");
	//ch_ideal->Add("ntuple_ideal_1E5.root");
	ch_ideal->Add("ntuple_ideal_2E5_seed1.root");
	ch_ideal->Add("ntuple_ideal_2E5_seed2.root");
	ch_ideal->Add("ntuple_ideal_2E5_seed3.root");

	// sanity check
	printf("full ntuple: %llu\n", ch_full->GetEntries());
	printf("ideal ntuple: %llu\n", ch_ideal->GetEntries());
	if (ch_full->GetEntries() != ch_ideal->GetEntries())
	{
		printf("ERROR: full and ideal ntuples are not compatible.\n");
		return 1;
	}
	
	// prepare output
	TFile *f_out = new TFile("analysis.root", "recreate");

	// select and link input branches
	ch_full->SetBranchStatus("*", 0);
	ch_ideal->SetBranchStatus("*", 0);

	EventMetaData *metaData = new EventMetaData();
	ch_full->SetBranchStatus("event_info.*", 1);
	ch_full->SetBranchAddress("event_info.", &metaData);

	//TriggerData *triggerData = new TriggerData();
	//ch_full->SetBranchStatus("trigger_data.*", 1);
	//ch_full->SetBranchAddress("trigger_data.", &triggerData);

	RPRootDumpReconstructedProton *simProton_R = new RPRootDumpReconstructedProton();
	ch_full->SetBranchStatus("sim_prot_right.*", 1);
	ch_full->SetBranchAddress("sim_prot_right.", &simProton_R);

	StationStruct st_L, st_R;
	st_L.AssignBranches("L", ch_full, ch_ideal, 20, 21, 24, 25);
	st_R.AssignBranches("R", ch_full, ch_ideal, 120, 121, 124, 125);

	// prepare counters and histograms
	map<string, CounterMap> counters;	// map: station label -> CounterMap
	
	// loop over events
	for (unsigned int en = 0; en < ch_full->GetEntries(); en++)
	{
		// load all data for enent en
		ch_full->GetEvent(en);
		ch_ideal->GetEvent(en);

		// get event number
		unsigned long ev = metaData->event_no;

		double th_y = simProton_R->thy;
		
		// run analysis
		st_L.Analyze(ev, th_y, counters["L"]);
		st_R.Analyze(ev, th_y, counters["R"]);
	}

	// print results
	for (map<string, CounterMap>::iterator stit = counters.begin(); stit != counters.end(); ++stit)
	{
		printf("\n\n==================== %s ====================\n", stit->first.c_str());

		for (map<string, map<string, unsigned long > >::iterator elit = stit->second.begin(); elit != stit->second.end(); ++elit)
		{
			printf("* %s\n", elit->first.c_str());

			for (map<string, unsigned long >::iterator cit = elit->second.begin(); cit != elit->second.end(); ++cit)
			{
				printf("\t%s: %lu\n", cit->first.c_str(), cit->second);
			}
		}

		CounterMap &c = stit->second;

		printf("\n");
		printf("N_T: N(prot & tr_val) / N(prot) = %.1f %%, N(prot & !tr_val) / N(prot) = %.1f %%\n",
			100. * c["N_T"]["prot & tr_val"] / c["N_T"]["prot"], 100. * c["N_T"]["prot & !tr_val"] / c["N_T"]["prot"]);
		printf("N_B: N(prot & tr_val) / N(prot) = %.1f %%, N(prot & !tr_val) / N(prot) = %.1f %%\n",
			100. * c["N_B"]["prot & tr_val"] / c["N_B"]["prot"], 100. * c["N_B"]["prot & !tr_val"] / c["N_B"]["prot"]);
		printf("F_T: N(prot & tr_val) / N(prot) = %.1f %%, N(prot & !tr_val) / N(prot) = %.1f %%\n",
			100. * c["F_T"]["prot & tr_val"] / c["F_T"]["prot"], 100. * c["F_T"]["prot & !tr_val"] / c["F_T"]["prot"]);
		printf("F_B: N(prot & tr_val) / N(prot) = %.1f %%, N(prot & !tr_val) / N(prot) = %.1f %%\n",
			100. * c["F_B"]["prot & tr_val"] / c["F_B"]["prot"], 100. * c["F_B"]["prot & !tr_val"] / c["F_B"]["prot"]);

		printf("\n");
		printf("N_T, F_T: N(prot & prot & tr_val & tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_T,F_T"]["prot & prot & tr_val & tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N_T, F_T: N(prot & prot & tr_val & !tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_T,F_T"]["prot & prot & tr_val & !tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N_T, F_T: N(prot & prot & !tr_val & tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_T,F_T"]["prot & prot & !tr_val & tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N_T, F_T: N(prot & prot & !tr_val & !tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_T,F_T"]["prot & prot & !tr_val & !tr_val"] / c["N_T,F_T"]["prot & prot"]);

		printf("\n");
		printf("N_B, F_B: N(prot & prot & tr_val & tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_B,F_B"]["prot & prot & tr_val & tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N_B, F_B: N(prot & prot & tr_val & !tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_B,F_B"]["prot & prot & tr_val & !tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N_B, F_B: N(prot & prot & !tr_val & tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_B,F_B"]["prot & prot & !tr_val & tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N_B, F_B: N(prot & prot & !tr_val & !tr_val) / N(prot & prot) = %.1f %%\n", 100. * c["N_B,F_B"]["prot & prot & !tr_val & !tr_val"] / c["N_B,F_B"]["prot & prot"]);
	
		printf("\n");
		printf("N(N_T.prot & F_T.prot & N_B.tr_val) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,N_B"]["prot & prot & opp-tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N(N_T.prot & F_T.prot & F_B.tr_val) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,F_B"]["prot & prot & opp-tr_val"] / c["N_T,F_T"]["prot & prot"]);
		printf("N(N_T.prot & F_T.prot & N_B.tr_val & F_B.tr_val) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,N_B,F_B"]["prot & prot & opp-tr_val & opp-tr_val"] / c["N_T,F_T"]["prot & prot"]);
	
		printf("\n");
		printf("N(N_T.prot & F_T.prot & N_B.pl_suff) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,N_B"]["prot & prot & opp-pl_suff"] / c["N_T,F_T"]["prot & prot"]);
		printf("N(N_T.prot & F_T.prot & F_B.pl_suff) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,F_B"]["prot & prot & opp-pl_suff"] / c["N_T,F_T"]["prot & prot"]);
		printf("N(N_T.prot & F_T.prot & N_B.pl_suff & F_B.pl_suff) / N(N_T.prot & F_T.prot) = %.3f %%\n",
			100. * c["N_T,F_T,N_B,F_B"]["prot & prot & opp-pl_suff & opp-pl_suff"] / c["N_T,F_T"]["prot & prot"]);
	
		printf("\n");
		printf("N(N_B.prot & F_B.prot & N_T.tr_val) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,N_T"]["prot & prot & opp-tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N(N_B.prot & F_B.prot & F_T.tr_val) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,F_T"]["prot & prot & opp-tr_val"] / c["N_B,F_B"]["prot & prot"]);
		printf("N(N_B.prot & F_B.prot & N_T.tr_val & F_T.tr_val) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,N_T,F_T"]["prot & prot & opp-tr_val & opp-tr_val"] / c["N_B,F_B"]["prot & prot"]);
	
		printf("\n");
		printf("N(N_B.prot & F_B.prot & N_T.pl_suff) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,N_T"]["prot & prot & opp-pl_suff"] / c["N_B,F_B"]["prot & prot"]);
		printf("N(N_B.prot & F_B.prot & F_T.pl_suff) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,F_T"]["prot & prot & opp-pl_suff"] / c["N_B,F_B"]["prot & prot"]);
		printf("N(N_B.prot & F_B.prot & N_T.pl_suff & F_T.pl_suff) / N(N_B.prot & F_B.prot) = %.3f %%\n",
			100. * c["N_B,F_B,N_T,F_T"]["prot & prot & opp-pl_suff & opp-pl_suff"] / c["N_B,F_B"]["prot & prot"]);
	}

	// save results
	st_L.MakeRatios();
	st_R.MakeRatios();

	st_L.Write();
	st_R.Write();

	delete f_out;
	return 0;
}
