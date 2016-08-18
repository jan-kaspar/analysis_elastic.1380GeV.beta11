#include "input_files.h"

#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"

#include <vector>
#include <map>
#include <string>

#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"
#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/RawDataFormats.h"
#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/RPRootTrackInfo.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

unsigned int verbosity;

//----------------------------------------------------------------------------------------------------

struct RPStruct
{
	RPRootDumpDigiInfo *digi;
	//RPRootDumpPatternInfo *pat;
	//RPRootDumpTrackInfo *tr;
	//vector<RPRootDumpTrackInfo> *mtr;

	TH1D *h_avg_mult;
	map<unsigned int, TH2D *> h2_avg_mult_fl_pl;

	RPStruct()
	{
		h_avg_mult = new TH1D("", ";avg mult", 51, -0.5, 50.5);

		for (unsigned int i = 2; i <= 4; i++)
			h2_avg_mult_fl_pl[i] = new TH2D("", ";avg mult first;avg mult last", 31, -0.5, 30.5, 31, -0.5, 30.5);
	}

	void AssignBranches(TChain *ch, unsigned int id)
	{
		char buf[100];

		digi = new RPRootDumpDigiInfo;
		sprintf(buf, "digi_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "digi_rp_%u.", id); ch->SetBranchAddress(buf, &digi);

		/*
		pat = new RPRootDumpPatternInfo();
		sprintf(buf, "nonpar_patterns_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "nonpar_patterns_rp_%u.", id); ch->SetBranchAddress(buf, &pat);

		tr = new RPRootDumpTrackInfo();
		sprintf(buf, "track_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch->SetBranchAddress(buf, &tr);

		mtr = new vector<RPRootDumpTrackInfo>();
		sprintf(buf, "multi_track_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "multi_track_rp_%u", id); ch->SetBranchAddress(buf, &mtr);
		*/
	}

	double Analyze();

	void Write(const char *name)
	{
		TDirectory *topDir = gDirectory;
		gDirectory = topDir->mkdir(name);
		
		h_avg_mult->Write("h_avg_mult");
	
		for (map<unsigned int, TH2D *>::iterator it = h2_avg_mult_fl_pl.begin(); it != h2_avg_mult_fl_pl.end(); ++it)
		{
			char buf[30];
			sprintf(buf, "h2_avg_mult_fl_pl_%u", it->first);
			it->second->Write(buf);
		}

		gDirectory = topDir;
	}
};

//----------------------------------------------------------------------------------------------------

double RPStruct::Analyze()
{
	if (verbosity > 1)
	{
		printf("\n\n>> RPStruct::Analyze\n");
	
		for (unsigned int i = 0; i < 10; i++)
			printf("\tplane %u: %u clusters\n", i, digi->numberOfClusters[i]);
	}

	// simple avergage
	double avg_mult = 0.;
	for (unsigned int i = 0; i < 10; i++)
		avg_mult += digi->numberOfClusters[i];
	avg_mult /= 10;

	if (verbosity > 1)
		printf("simple avg_mult = %f\n", avg_mult);

	// sort planes according to multiplicity difference from the average
	vector< pair<double, unsigned int> > list;
	for (unsigned int i = 0; i < 10; i++)
		list.push_back({fabs(digi->numberOfClusters[i] - avg_mult), i});
	sort(list.begin(), list.end(),
			[](const pair<double, int> &left, const pair<double, int> &right)
			{
    			return left.first < right.first;
			}
		);

	// recalculte average over a given number of planes with the lowest difference
	unsigned int N_keep_planes = 8;
	avg_mult = 0.;
	for (unsigned int i = 0; i < N_keep_planes; i++)
		avg_mult += digi->numberOfClusters[i];
	avg_mult /= N_keep_planes;

	if (verbosity > 1)
		printf("refined avg_mult = %f\n", avg_mult);

	// fill average multiplicity
	h_avg_mult->Fill(avg_mult);

	// average multiplicity in first and last planes
	for (map<unsigned int, TH2D *>::iterator it = h2_avg_mult_fl_pl.begin(); it != h2_avg_mult_fl_pl.end(); ++it)
	{
		unsigned int N_pl = it->first;
		double S_f = 0., S_l = 0.;
		for (unsigned int i = 0; i < N_pl; i++)
		{
			S_f += digi->numberOfClusters[i];
			S_l += digi->numberOfClusters[9 - i];
		}

		it->second->Fill(S_f/N_pl, S_l/N_pl);
		
		if (verbosity > 1)
			printf("\t%u-plane averages: first=%.3f, last=%.3f\n", N_pl, S_f/N_pl, S_l/N_pl);
	}

	return avg_mult;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct StationStruct
{
	RPStruct N_T, N_B, F_T, F_B;

	vector< vector<unsigned int> > configurations;
	vector< double > thresholds;
	vector< pair<TH1D*, TH1D*> > histograms;

	StationStruct()
	{
		configurations.push_back({0, 4});
		configurations.push_back({1, 5});
		configurations.push_back({0, 4, 5});
		configurations.push_back({1, 4, 5});
		configurations.push_back({0, 1, 4, 5});
		configurations.push_back({0, 1, 4});
		configurations.push_back({0, 1, 5});

		for (unsigned int i = 1; i <= 20; i++)
			thresholds.push_back(i);

		for (unsigned int i = 0; i < configurations.size(); i++)
		{
			TH1D *h_inc = new TH1D("", ";threshold", 21, -0.5, 20.5);
			TH1D *h_exc = new TH1D("", ";threshold", 21, -0.5, 20.5);
			histograms.push_back({h_inc, h_exc});
		}
	}

	void AssignBranches(TChain *ch, unsigned int nt, unsigned int nb, unsigned int ft, unsigned int fb)
	{
		N_T.AssignBranches(ch, nt);
		N_B.AssignBranches(ch, nb);
		F_T.AssignBranches(ch, ft);
		F_B.AssignBranches(ch, fb);
	}

	void Analyze(double, double, double, double);

	void Write(const char *name)
	{
		TDirectory *topDir = gDirectory;
		TDirectory *stDir = topDir->mkdir(name);

		gDirectory = stDir;
		N_T.Write("N_T");
		N_B.Write("N_B");
		F_T.Write("F_T");
		F_B.Write("F_B");
		
		for (unsigned int i = 0; i < configurations.size(); i++)
		{
			string label;
			for (unsigned int rpi = 0; rpi < configurations[i].size(); rpi++)
			{
				if (rpi > 0) label += ",";
				if (configurations[i][rpi] == 0) label += "N_T";
				if (configurations[i][rpi] == 1) label += "N_B";
				if (configurations[i][rpi] == 4) label += "F_T";
				if (configurations[i][rpi] == 5) label += "F_B";
			}

			gDirectory = stDir->mkdir(label.c_str());
			histograms[i].first->Write("h_inc");
			histograms[i].second->Write("h_exc");
		}

		gDirectory = topDir;
	}
};

//----------------------------------------------------------------------------------------------------

void StationStruct::Analyze(double amNT, double amNB, double amFT, double amFB)
{
	if (verbosity > 1)
		printf("\n\n>> StationStruct::Analyze(NT=%.2f, NB=%.2f, FT=%.2f, FB=%.2f)\n", amNT, amNB, amFT, amFB);

	map<unsigned int, double> am_map;
	am_map[0] = amNT;
	am_map[1] = amNB;
	am_map[4] = amFT;
	am_map[5] = amFB;

	for (unsigned int ci = 0; ci < configurations.size(); ci++)
	{
		for (unsigned int ti = 0; ti < thresholds.size(); ti++)	
		{
			double th = thresholds[ti];

			bool all_sel_above = true, all_non_sel_below = true;

			for (map<unsigned int, double>::iterator it = am_map.begin(); it != am_map.end(); ++it)
			{
				bool sel = (find(configurations[ci].begin(), configurations[ci].end(), it->first) != configurations[ci].end());

				if (sel)
					all_sel_above &= (it->second >= th);
				else
					all_non_sel_below &= (it->second < th);
			}

			if (all_sel_above)
			{
				histograms[ci].first->Fill(th);
				
				if (verbosity > 1)
					printf("\tfill: configuration %u, threshold %.2f, inc\n", ci, th);
			}

			if (all_sel_above && all_non_sel_below)
			{
				histograms[ci].second->Fill(th);

				if (verbosity > 1)
					printf("\tfill: configuration %u, threshold %.2f, exc\n", ci, th);
			}
		}
	}
}

//----------------------------------------------------------------------------------------------------

void AnalyzeStation(StationStruct &st)
{
	unsigned int mult_N_T = st.N_T.Analyze();
	unsigned int mult_N_B = st.N_B.Analyze();
	unsigned int mult_F_T = st.F_T.Analyze();
	unsigned int mult_F_B = st.F_B.Analyze();

	st.Analyze(mult_N_T, mult_N_B, mult_F_T, mult_F_B);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// configuration
	verbosity = 0;

	// init diagonal
	Init(argv[1]);
	if (diagonal != dCombined)
		return rcIncompatibleDiagonal;

	// get input
	InitInputFiles();
	TChain *ch = new TChain("TotemNtuple");
	for (unsigned int i = 0; i < input_files.size(); i++)
	{
		ch->Add(input_files[i].c_str());
		printf("+ %s\n", input_files[i].c_str());
	}
	
	// prepare output
	TFile *outF = new TFile((string("showers_") + argv[1] + ".root").c_str(), "recreate");

	// select and link input branches
	ch->SetBranchStatus("*", 0);

	EventMetaData *metaData = new EventMetaData();
	ch->SetBranchStatus("event_info.*", 1);
	ch->SetBranchAddress("event_info.", &metaData);

	TriggerData *triggerData = new TriggerData();
	ch->SetBranchStatus("trigger_data.*", 1);
	ch->SetBranchAddress("trigger_data.", &triggerData);

	StationStruct st_L, st_R;
	st_L.AssignBranches(ch, 20, 21, 24, 25);
	st_R.AssignBranches(ch, 120, 121, 124, 125);

	signed int offset = 0;
	unsigned int prev_run = 0;

	// loop over events
	for (unsigned int ev = 0; ev < ch->GetEntries(); ev++)
	{
		// TODO: remove
		//if (ev > 1000)
		//	break;

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
		
		// zero-bias selection
		if (! IsZeroBias(td.input_status_bits, metaData->run_no, metaData->event_no))
			continue;

		// as zero-bias data are used, blinding (some) elastic events is not needed

		// skip troublesome runs
		if (SkipRun(run, file, false))
			continue;

		// select bunches
		if (SkipBunch(run, td.bunch_num))
			continue;

		// run analysis
		AnalyzeStation(st_L);
		AnalyzeStation(st_R);
	}

	// save results
	st_L.Write("L");
	st_R.Write("R");

	delete outF;
	return 0;
}
