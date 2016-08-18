#include "TChain.h"

#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"
#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/RawDataFormats.h"
#include "/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAnalysis/TotemNtuplizer/interface/RPRootTrackInfo.h"

#include <set>
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct Event
{
	unsigned int num;
	string desc;
};

bool operator<(const Event &l, const Event &r)
{
	return l.num < r.num;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// settings
	set<Event> event_selection;

        event_selection.insert({129, "L, N_B & F_B"});
        event_selection.insert({266, "R, N_B & F_B"});
        event_selection.insert({370, "L, N_B & F_B"});
        event_selection.insert({546, "L, N_B & F_B"});
        event_selection.insert({616, "L, N_B & F_B"});
        event_selection.insert({865, "L, N_B & F_B"});
        event_selection.insert({1055, "L, N_B & F_B"});
        event_selection.insert({1243, "R, N_T & F_T"});
        event_selection.insert({1318, "L, N_T & F_T"});
        event_selection.insert({1419, "L, N_T & F_T"});
        event_selection.insert({1478, "L, N_B & F_B"});
        event_selection.insert({1637, "L, N_T & F_T"});
        event_selection.insert({1778, "R, N_B & F_B"});
        event_selection.insert({1804, "R, N_B & F_B"});
        event_selection.insert({1846, "R, N_B & F_B"});
        event_selection.insert({1882, "R, N_T & F_T"});
        event_selection.insert({1906, "R, N_B & F_B"});
        event_selection.insert({2240, "R, N_B & F_B"});
        event_selection.insert({2365, "R, N_B & F_B"});
        event_selection.insert({2393, "R, N_B & F_B"});
        event_selection.insert({2420, "R, N_B & F_B"});
        event_selection.insert({2423, "L, N_T & F_T"});

	// get input
	TChain *ch = new TChain("ntuple");
	ch->Add("ntuple_sim_hits.root");

	TChain *ch_full = new TChain("TotemNtuple");
	ch_full->Add("ntuple_full.root");

	TChain *ch_ideal = new TChain("TotemNtuple");
	ch_ideal->Add("ntuple_ideal.root");

	//printf("entries: %llu\n", ch->GetEntries());
	
	// assign branches - sim hits
	ch->SetBranchStatus("*", 0);

	double event, energyLoss, trackId, particleType;
	double x, y, z;

	ch->SetBranchStatus("Event", 1); ch->SetBranchAddress("Event", &event);
	ch->SetBranchStatus("x", 1); ch->SetBranchAddress("x", &x);
	ch->SetBranchStatus("y", 1); ch->SetBranchAddress("y", &y);
	ch->SetBranchStatus("z", 1); ch->SetBranchAddress("z", &z);
	ch->SetBranchStatus("ELoss", 1); ch->SetBranchAddress("ELoss", &energyLoss);
	ch->SetBranchStatus("TrackID", 1); ch->SetBranchAddress("TrackID", &trackId);
	ch->SetBranchStatus("Ptype", 1); ch->SetBranchAddress("Ptype", &particleType);

	// assign branches - RP tracks
	ch_full->SetBranchStatus("*", 0);
	ch_ideal->SetBranchStatus("*", 0);

	vector<unsigned int> RPs = { 20, 21, 24, 25, 120, 121, 124, 125 };

	map<unsigned int, RPRootDumpTrackInfo*> rp_track, rp_track_id;
	map<unsigned int, RPRootDumpDigiInfo*> rp_digi, rp_digi_id;

	for (unsigned int i = 0; i < RPs.size(); i++)
	{
		rp_track[RPs[i]] = new RPRootDumpTrackInfo();
		rp_track_id[RPs[i]] = new RPRootDumpTrackInfo();

		rp_digi[RPs[i]] = new RPRootDumpDigiInfo();
		rp_digi_id[RPs[i]] = new RPRootDumpDigiInfo();
	}

	for (unsigned int i = 0; i < RPs.size(); i++)
	{
		char buf[100];
		unsigned int id = RPs[i];

		sprintf(buf, "track_rp_%u.*", id); ch_full->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch_full->SetBranchAddress(buf, &rp_track[id]);

		sprintf(buf, "track_rp_%u.*", id); ch_ideal->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch_ideal->SetBranchAddress(buf, &rp_track_id[id]);

		sprintf(buf, "digi_rp_%u.*", id); ch_full->SetBranchStatus(buf, 1);
		sprintf(buf, "digi_rp_%u.", id); ch_full->SetBranchAddress(buf, &rp_digi[id]);

		//sprintf(buf, "digi_rp_%u.*", id); ch_ideal->SetBranchStatus(buf, 1);
		//sprintf(buf, "digi_rp_%u.", id); ch_ideal->SetBranchAddress(buf, &rp_digi_id[id]);
	}

	// process input
	bool prev_ev_set = false;
	unsigned int prev_ev = 0;
	unsigned int entries = ch->GetEntries();
	bool active = false;
	for (unsigned int en = 0; en < entries; en++)
	{
		ch->GetEvent(en);

		unsigned int ev = (unsigned int) floor(event + 0.5);

		// delimit event data
		bool event_change = !prev_ev_set || (prev_ev != ev);

		if (event_change && active)
		{
			// dump RP reco data
			int reco_entry = prev_ev - 1;
			ch_full->GetEvent(reco_entry);
			ch_ideal->GetEvent(reco_entry);

			for (unsigned int rpi = 0; rpi < RPs.size(); rpi++)
			{
				unsigned int rpId = RPs[rpi];

				RPRootDumpTrackInfo *t = rp_track[rpId];
				if (t->valid)
					printf("TrackFit(%u, %e, %e, %e, %e, %e, true);\n", rpId, t->x, t->y, t->z, t->thx, t->thy);

				t = rp_track_id[rpId];
				if (t->valid)
					printf("TrackFit(%u, %e, %e, %e, %e, %e, false);\n", rpId, t->x, t->y, t->z, t->thx, t->thy);

				RPRootDumpDigiInfo *d = NULL;

				d = rp_digi[rpId];
				if (d->numberOfPlanesOn > 0)
				{
					printf("RPClusters(%u, true", rpId);
					for (unsigned int i = 0; i < 10; i++)
						printf(", %i", d->numberOfClusters[i]);
					printf(");\n");
				}
			}


			printf("EndEvent(%u);\n", prev_ev);
			active = false;
		}

		prev_ev = ev;
		prev_ev_set = true;

		// event selected for dump?
		set<Event>::iterator evIt = event_selection.find({ev, ""});
		bool ev_selected = (evIt != event_selection.end());

		if (ev_selected)
		{
			if (event_change)
			{
				printf("\nBeginEvent(%u, \"%s\");\n", ev, evIt->desc.c_str());
				active = true;
			}

			bool interesting_hit = false;
			if (fabs(z) > 214.578E3 && fabs(z) < 214.678E3) interesting_hit = true;
			if (fabs(z) > 219.950E3 && fabs(z) < 220.050E3) interesting_hit = true;

			// print event data
			if (interesting_hit)
				printf("Hit(%e, %e, %e, %e, %.0f, %.0f);\n", x, y, z, energyLoss, trackId, particleType);
		}
	}

	return 0;
}
