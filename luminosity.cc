#include "../common_definitions.h"
#include "../common_algorithms.h"
#include "parameters.h"
#include "../common.h"
//#include "../LumiCalculator.h"

#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TProfile.h"

#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct IndexTimePair
{
	unsigned int index, timestamp;
	IndexTimePair(unsigned int _i, unsigned int _t) : index(_i), timestamp(_t) {}

	bool operator< (const IndexTimePair &other) const
	{
		return (timestamp < other.timestamp);
	}
};

//----------------------------------------------------------------------------------------------------

struct period
{
	double first, last;
	unsigned int ev_first, ev_last;
	unsigned int tr_first, tr_last;
	unsigned int ev_total;
	unsigned int run;

	double daq_eff;
	double L_int;

	period() : first(0.), last(0.), ev_first(0), ev_last(0), tr_first(0), tr_last(0), ev_total(0), run(0), daq_eff(0.), L_int(0.) {}
};

//----------------------------------------------------------------------------------------------------

// bunch number --> graph: time vs luminosity
map<unsigned int, TGraph*> lumiData;

void LoadLumiData(const string &file, const vector<unsigned int> bunches)
{
	FILE *f = fopen(file.c_str(), "r");
	if (!f)
	{
		printf("ERROR in LoadLumiData > Can't read file `%s'.\n", file.c_str());
		return;
	}

	while(!feof(f))
	{
		char buf[100000];
		buf[0] = 0;
		fgets(buf, 100000, f);

		if (strlen(buf) == 0)
			continue;
		
		if (buf[0] == '*' || buf[0] == '#')
			continue;

		char *p;

		p = strtok(buf, ",");
		if (p == NULL)
			break;

		string runFill = p;

		p = strtok(NULL, ",");
		string lumiSec = p;
		
		p = strtok(NULL, ",");
		string time = p;

		tm time_tm;
		strptime(time.c_str(), "%m/%d/%y %H:%M:%S", &time_tm);
		time_t timestamp = mktime(&time_tm) - timestamp0;
		//printf("%s => %u\n", time.c_str(), timestamp);

		// TODO: applicable ???
		// correction: UTC --> Central European Summer Time
		//timestamp += 2*3600;
		
		p = strtok(NULL, ",");
		string fullDelivered = p;
		
		p = strtok(NULL, ",");
		string fullRecorded = p;

		while (true)
		{
			p = strtok(NULL, ",");
			if (!p)
				break;
			unsigned int bunch = atoi(p);

			p = strtok(NULL, ",");
			if (!p)
				break;
			double lumi = atof(p);
			
			if (find(bunches.begin(), bunches.end(), bunch) != bunches.end())
			{
				TGraph *g = lumiData[bunch];
				if (!g)
				{
					g = lumiData[bunch] = new TGraph();
					char name[20];
					sprintf(name,"lumi_%u", bunch);
					g->SetName(name);
				}

				int idx = g->GetN();
				g->SetPoint(idx, timestamp, lumi);
			}
		}
	}

	fclose(f);
}

//----------------------------------------------------------------------------------------------------

void SaveLumiData()
{
	for (map<unsigned int, TGraph*>::iterator it = lumiData.begin(); it != lumiData.end(); ++it)
	{
		it->second->Write();
	}
}

//----------------------------------------------------------------------------------------------------

double CalculateLumi(unsigned int run, unsigned int first, unsigned int last)
{
	// determine active bunches (TOTEM numbering scheme)
	const vector<unsigned int> &bunches = bunchMap[run];

	if (bunches.size() == 0)
	{
		printf("ERROR: no bunches for run %u.\n", run);
	}

	double S_all_bunch = 0.;
	for (unsigned int bi = 0; bi < bunches.size(); ++bi)
	{
		// bunch in LHC numbering scheme
		unsigned int bunch_num_lhc = bunches[bi] + 1;

		// luminosity data
		TGraph *g_lumi = lumiData[bunch_num_lhc];
		if (! g_lumi)
		{
			printf("ERROR: g_lumi = NULL\n");
			continue;
		}
	
		// IMPORTANT: assuming negligible time skew (verified for runs 8333 to 8372)
	
		// collect integration data
		vector<pair<unsigned int, double> > data;
		data.push_back(pair<unsigned int, double>(first, g_lumi->Eval(first)));
		
		for (int i = 0; i < g_lumi->GetN(); i++)
		{
			double time, l;
			g_lumi->GetPoint(i, time, l);
			unsigned int time_i = floor(time + 0.5);
			if (time_i > first && time_i < last)
				data.push_back(pair<unsigned int, double>(time, l));
		}
		
		data.push_back(pair<unsigned int, double>(last, g_lumi->Eval(last)));
	
		/*
		printf("bunch %u: \n", bunch_num_lhc);
		printf("int. points %lu\n", data.size());
	
		printf("first = %u\n", first);
		printf("last = %u\n", last);
	
		for (unsigned int i = 0; i < data.size(); i++)
			printf("\ttime = %u, lumi = %f\n", data[i].first, data[i].second);
		*/
	
		// do integration
		double S = 0.;
		for (unsigned int i = 1; i < data.size(); i++)
		{
			S += (data[i].second + data[i-1].second) / 2. * (data[i].first - data[i-1].first);
		}
	
		//printf("mean lumi %E\n", S / (last - first));

		if (bi > 0)
			printf(" + ");
		printf("%.3E (%u)", S * 1E3, bunch_num_lhc);	// conversion to mb^-1

		S_all_bunch += S;
	}

	return S_all_bunch;
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

	// load lumi data
	vector<unsigned int> selected_bunches; // our bunch 0 corresponds to bunch 1 in the CMS notation
	selected_bunches.push_back(0); 
	selected_bunches.push_back(1); 

	printf(">> luminosity_data_file: %s\n", luminosity_data_file.c_str());
	LoadLumiData(luminosity_data_file, selected_bunches);

	// print conventions
	printf(">> timestamp = UNIX time - %.0f\n", timestamp0);
	
	// get input data
	TFile *distF = new TFile((string("distributions_") + argv[1] + ".root").c_str());
	TGraph *g_timestamp_vs_ev_idx_sel = (TGraph *) distF->Get("metadata/g_timestamp_vs_ev_idx_sel");

	TFile *dataF = new TFile((string("distill_") + argv[1] + ".root").c_str());
	TTree *inT = (TTree *) dataF->Get("distilled");
	EventRed ev;
	
	inT->SetBranchAddress("timestamp", &ev.timestamp);
	inT->SetBranchAddress("run_num", &ev.run_num);
	inT->SetBranchAddress("event_num", &ev.event_num);
	inT->SetBranchAddress("trigger_num", &ev.trigger_num);

	// prepare output
	TFile *outF = new TFile((string("luminosity_") + argv[1] + ".root").c_str(), "recreate");
	SaveLumiData();

	// prepare plots
	TH1D *h_daqEff = new TH1D("h_daqEff", ";DAQ efficiency", 50, 0.95, 1.00); h_daqEff->SetLineColor(2);

	// sort list of selected indeces according to time
	vector<IndexTimePair> indeces;
	for (int si = 0; si < g_timestamp_vs_ev_idx_sel->GetN(); si++)
	{
		double ri, rt;
		g_timestamp_vs_ev_idx_sel->GetPoint(si, ri, rt);
		indeces.push_back(IndexTimePair(ri, rt));
	}

	sort(indeces.begin(), indeces.end());

	// find periods
	vector<period> periods;
	for (unsigned long si = 0; si < indeces.size(); si++)
	{
		// use the same event selection as in distributions.root
		unsigned int ev_idx = indeces[si].index;
		inT->GetEntry(ev_idx);

		// TODO: remove debug
		/*
		//printf("%u, %u\n", indeces[si].timestamp, indeces[si].index);
		if (si > 200000)
			break;
		*/

		unsigned int run = ev.run_num / 10000;
		//unsigned int run = ev.run_num;

		// find an existing and compatible period
		bool found = false;
		for (unsigned int pi = 0; pi < periods.size(); pi++)
		{
			unsigned int margin = 3;
			if (periods[pi].run == run && ev.timestamp >= periods[pi].first - margin && ev.timestamp <= periods[pi].last + margin)
			{
				found = true;
				//printf("\tFOUND %u\n", pi);

				periods[pi].first = min(periods[pi].first, double(ev.timestamp));
				periods[pi].last = max(periods[pi].last, double(ev.timestamp));
				
				periods[pi].ev_first = min(periods[pi].ev_first, ev.event_num);
				periods[pi].ev_last = max(periods[pi].ev_last, ev.event_num);
				
				periods[pi].tr_first = min(periods[pi].tr_first, ev.trigger_num);
				periods[pi].tr_last = max(periods[pi].tr_last, ev.trigger_num);

				periods[pi].ev_total++;

				break;
			}
		}

		// add new period
		if (!found)
		{
			period p;
			p.first = ev.timestamp;
			p.last = ev.timestamp;
			
			p.ev_first = ev.event_num;
			p.ev_last = ev.event_num;
			
			p.tr_first = ev.trigger_num;
			p.tr_last = ev.trigger_num;

			p.ev_total = 1;
			p.run = run;

			periods.push_back(p);
		}
	}

	// prune unreasonable periods
	vector<period> selected_periods;
	for (unsigned int pi = 0; pi < periods.size(); pi++)
	{
		period &p = periods[pi];
		double duration = p.last - p.first;
		unsigned int events = p.ev_last - p.ev_first;

		bool keep = true;

		if (duration <= 10.)
		{
			printf("period %u: too little duration %.1f\n", pi, duration);
			keep = false;
		}
		
		if (events < 10)
		{
			printf("period %u: too little events %u\n", pi, events+1);
			keep = false;
		}
		
		if (keep)	
			selected_periods.push_back(p);
		else
			printf("\t=> removing period %u (from %.1f to %.1f, run %u)\n", pi, p.first, p.last, p.run);
	}

	// print and sum
	printf("\n");
	printf("per| run  | timestamp first, last : diff  | event cnt first, last : diff  | trigger cnt first, last : diff|sel. evnts| DAQ eff  |L_int (mb^-1)\n");
	
	double effLumiSum = 0.;
	double durationSum = 0.;
	double effWeightedSum = 0.;
	
	for (unsigned int pi = 0; pi < selected_periods.size(); pi++)
	{
		period &p = selected_periods[pi];

		double duration = p.last - p.first;
		unsigned int events = p.ev_last - p.ev_first;
		unsigned int triggers = p.tr_last - p.tr_first;

		p.daq_eff = double(events) / triggers;

		printf("%3u| %u | %8.1f, %8.1f : %8.1f | %8u, %8u : %8u | %8u, %8u : %8u | %8u | %8.3f | ",
			pi, p.run,
			p.first, p.last, duration,
			p.ev_first, p.ev_last, events,
			p.tr_first, p.tr_last, triggers,
			p.ev_total, p.daq_eff);
		
		p.L_int = CalculateLumi(p.run, p.first, p.last) * 1E3;	// conversion to mb^-1

		printf(" = %11.3E\n", p.L_int);
		
		effLumiSum += p.L_int * p.daq_eff;
		durationSum += duration;
		effWeightedSum += p.daq_eff * duration;

		h_daqEff->Fill(p.daq_eff);
	}

	printf("\n");
	printf(">> total duration: %.0f s\n", durationSum);
	printf(">> mean DAQ effeciency: %.3f\n", effWeightedSum / durationSum);
	printf(">> effective luminosity sum: %.3E mb^-1\n", effLumiSum);

	h_daqEff->Write();

	delete outF;

	return 0;
}
