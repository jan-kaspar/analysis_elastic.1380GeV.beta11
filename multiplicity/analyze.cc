#include "TFile.h"
#include "TChain.h"

#include <vector>
#include <map>
#include <string>

#include "TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RawDataFormats.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RPRootTrackInfo.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

struct RPStruct
{
	RPRootDumpDigiInfo *digi;
	/*
	RPRootDumpPatternInfo *pat;
	RPRootDumpTrackInfo *tr;
	*/
	vector<RPRootDumpTrackInfo> *mtr;

	RPStruct() : digi(NULL) /*, pat(NULL), tr(NULL) */ {}

	void AssignBranches(TChain *ch, unsigned int id)
	{
		char buf[100];

		digi = new RPRootDumpDigiInfo;
		sprintf(buf, "digi_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "digi_rp_%u.", id); ch->SetBranchAddress(buf, &digi);

		/*
		pat = new RPRootDumpPatternInfo();
		sprintf(buf, "patterns_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "patterns_rp_%u.", id); ch->SetBranchAddress(buf, &pat);

		tr = new RPRootDumpTrackInfo();
		sprintf(buf, "track_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "track_rp_%u.", id); ch->SetBranchAddress(buf, &tr);
		*/

		mtr = new vector<RPRootDumpTrackInfo>();
		sprintf(buf, "multi_track_rp_%u.*", id); ch->SetBranchStatus(buf, 1);
		sprintf(buf, "multi_track_rp_%u", id); ch->SetBranchAddress(buf, &mtr);
	}
};

//----------------------------------------------------------------------------------------------------

struct DiagStruct
{
	RPStruct L_2_F, L_2_N, R_2_N, R_2_F;

	void AssignBranches(TChain *ch, unsigned int l2f, unsigned int l2n, unsigned int r2n, unsigned int r2f)
	{
		L_2_F.AssignBranches(ch, l2f);
		L_2_N.AssignBranches(ch, l2n);
		R_2_N.AssignBranches(ch, r2n);
		R_2_F.AssignBranches(ch, r2f);
	}
};

//----------------------------------------------------------------------------------------------------

// map: label -> number of events
typedef map<string, unsigned int> CounterMap;

//----------------------------------------------------------------------------------------------------

enum Projection { pU, pV };

//----------------------------------------------------------------------------------------------------

bool ConditionFulfilled(const string &cond, unsigned int val)
{
	const char &op = cond.at(0);
	unsigned int cval = atoi(cond.substr(1).c_str());

	if (op == '=' && val == cval)
		return true;

	if (op == '<' && val < cval)
		return true;

	if (op == '>' && val > cval)
		return true;

	return false;
}

//----------------------------------------------------------------------------------------------------

// USAGE WARNING: put more restrictive conditions first, more lose at the end

bool FitsMultiplicityCondition(const vector<int> &multiplicities, Projection proj, vector<string> conditions)
{
	if (conditions.size() != 5)
		throw "condition size is not 5";

	//printf("------------------------------------------------\n");

	vector<unsigned int> plIndeces;
	if (proj == pU)
		plIndeces = { 1, 3, 5, 7, 9 };
   	if (proj == pV)
		plIndeces = { 0, 2, 4, 6, 8 };

	/*
	printf("multiplicities: ");
	for (const auto plIdx : plIndeces)
		printf("%u --> %u, ", plIdx, multiplicities[plIdx]);
	printf("\n");
	*/

	for (unsigned int ci = 0; ci < conditions.size(); ++ci)
	{
		//printf("ci = %u (%s)\n", ci, conditions[ci].c_str());

		bool planeFound = false;
		for (auto plIt = plIndeces.begin(); plIt != plIndeces.end(); ++plIt)
		{
			//printf("    plIt = %u\n", *plIt);

			if (ConditionFulfilled(conditions[ci], multiplicities[*plIt]))
			{
				//printf("        found: %u\n", multiplicities[*plIt]);
				plIndeces.erase(plIt);
				planeFound = true;
				break;
			}
		}

		if (!planeFound)
		{
			//printf("    no matching plane found\n");
			return false;
		}
	}

	return true;
}

//----------------------------------------------------------------------------------------------------

bool FitsMultiplicityCondition(const RPStruct &rp, Projection proj, vector<string> conditions)
{
	const auto &multiplicities = rp.digi->numberOfClusters;
	return FitsMultiplicityCondition(multiplicities, proj, conditions);
}

//----------------------------------------------------------------------------------------------------

bool IsReconstructable(const RPStruct &rp)
{
	unsigned int n_reas_V = 0;
	for (unsigned int idx = 0; idx < 10; idx += 2)
	{
		if (rp.digi->numberOfClusters[idx] <= 5)
			n_reas_V++;
	}

	unsigned int n_reas_U = 0;
	for (unsigned int idx = 1; idx < 10; idx += 2)
	{
		if (rp.digi->numberOfClusters[idx] <= 5)
			n_reas_U++;
	}

	if (n_reas_V >= 3 && n_reas_U >= 3)
		return true;

	return false;
}

//----------------------------------------------------------------------------------------------------

bool HasMultitrack(const RPStruct &rp)
{
	for (const auto &t : *rp.mtr)
	{
		if (t.valid)
			return true;
	}

	return false;
}

//----------------------------------------------------------------------------------------------------

void AnalyzeDiagonal(const DiagStruct &dgn, unsigned int /*evIdx*/, CounterMap &c)
{
	c["anything"]++;

	//--------------------

	bool L_2_F_mtr = HasMultitrack(dgn.L_2_F);
	bool L_2_N_mtr = HasMultitrack(dgn.L_2_N);
	bool R_2_N_mtr = HasMultitrack(dgn.R_2_N);
	bool R_2_F_mtr = HasMultitrack(dgn.R_2_F);

	if (L_2_F_mtr) c["L_2_F_mtr"]++;
	if (L_2_N_mtr) c["L_2_N_mtr"]++;
	if (R_2_N_mtr) c["R_2_N_mtr"]++;
	if (R_2_F_mtr) c["R_2_F_mtr"]++;

	bool all_rp_mtr = L_2_F_mtr && L_2_N_mtr && R_2_N_mtr && R_2_F_mtr;

	if (all_rp_mtr) c["all_rp_mtr"]++;

	// TODO
	/*
	if (!all_rp_mtr)
		return;
	*/

	//--------------------

	bool L_2_F_U_cond0 = FitsMultiplicityCondition(dgn.L_2_F, pU, {"=5", "=5", "=5", ">5", ">5"});
	bool L_2_F_V_cond0 = FitsMultiplicityCondition(dgn.L_2_F, pV, {"=5", "=5", "=5", ">5", ">5"});

	bool L_2_N_U_cond0 = FitsMultiplicityCondition(dgn.L_2_N, pU, {"=5", "=5", "=5", ">5", ">5"});
	bool L_2_N_V_cond0 = FitsMultiplicityCondition(dgn.L_2_N, pV, {"=5", "=5", "=5", ">5", ">5"});

	bool R_2_N_U_cond0 = FitsMultiplicityCondition(dgn.R_2_N, pU, {"=5", "=5", "=5", ">5", ">5"});
	bool R_2_N_V_cond0 = FitsMultiplicityCondition(dgn.R_2_N, pV, {"=5", "=5", "=5", ">5", ">5"});

	bool R_2_F_U_cond0 = FitsMultiplicityCondition(dgn.R_2_F, pU, {"=5", "=5", "=5", ">5", ">5"});
	bool R_2_F_V_cond0 = FitsMultiplicityCondition(dgn.R_2_F, pV, {"=5", "=5", "=5", ">5", ">5"});

	//--------------------

	bool L_2_F_U_cond1m = FitsMultiplicityCondition(dgn.L_2_F, pU, {"=4", "=5", "=5", ">5", ">5"});
	bool L_2_F_V_cond1m = FitsMultiplicityCondition(dgn.L_2_F, pV, {"=4", "=5", "=5", ">5", ">5"});

	bool L_2_N_U_cond1m = FitsMultiplicityCondition(dgn.L_2_N, pU, {"=4", "=5", "=5", ">5", ">5"});
	bool L_2_N_V_cond1m = FitsMultiplicityCondition(dgn.L_2_N, pV, {"=4", "=5", "=5", ">5", ">5"});

	bool R_2_N_U_cond1m = FitsMultiplicityCondition(dgn.R_2_N, pU, {"=4", "=5", "=5", ">5", ">5"});
	bool R_2_N_V_cond1m = FitsMultiplicityCondition(dgn.R_2_N, pV, {"=4", "=5", "=5", ">5", ">5"});

	bool R_2_F_U_cond1m = FitsMultiplicityCondition(dgn.R_2_F, pU, {"=4", "=5", "=5", ">5", ">5"});
	bool R_2_F_V_cond1m = FitsMultiplicityCondition(dgn.R_2_F, pV, {"=4", "=5", "=5", ">5", ">5"});

	//--------------------

	bool L_2_F_U_cond1p = FitsMultiplicityCondition(dgn.L_2_F, pU, {"=5", "=5", "=6", ">5", ">5"});
	bool L_2_F_V_cond1p = FitsMultiplicityCondition(dgn.L_2_F, pV, {"=5", "=5", "=6", ">5", ">5"});

	bool L_2_N_U_cond1p = FitsMultiplicityCondition(dgn.L_2_N, pU, {"=5", "=5", "=6", ">5", ">5"});
	bool L_2_N_V_cond1p = FitsMultiplicityCondition(dgn.L_2_N, pV, {"=5", "=5", "=6", ">5", ">5"});

	bool R_2_N_U_cond1p = FitsMultiplicityCondition(dgn.R_2_N, pU, {"=5", "=5", "=6", ">5", ">5"});
	bool R_2_N_V_cond1p = FitsMultiplicityCondition(dgn.R_2_N, pV, {"=5", "=5", "=6", ">5", ">5"});

	bool R_2_F_U_cond1p = FitsMultiplicityCondition(dgn.R_2_F, pU, {"=5", "=5", "=6", ">5", ">5"});
	bool R_2_F_V_cond1p = FitsMultiplicityCondition(dgn.R_2_F, pV, {"=5", "=5", "=6", ">5", ">5"});

	//--------------------
	
	bool L_2_F_recon = IsReconstructable(dgn.L_2_F);
	bool L_2_N_recon = IsReconstructable(dgn.L_2_N);
	bool R_2_N_recon = IsReconstructable(dgn.R_2_N);
	bool R_2_F_recon = IsReconstructable(dgn.R_2_F);

	//--------------------

	if (L_2_F_U_cond0) c["L_2_F_U_cond0"]++;
	if (L_2_F_V_cond0) c["L_2_F_V_cond0"]++;

	if (L_2_N_U_cond0) c["L_2_N_U_cond0"]++;
	if (L_2_N_V_cond0) c["L_2_N_V_cond0"]++;

	if (R_2_N_U_cond0) c["R_2_N_U_cond0"]++;
	if (R_2_N_V_cond0) c["R_2_N_V_cond0"]++;

	if (R_2_F_U_cond0) c["R_2_F_U_cond0"]++;
	if (R_2_F_V_cond0) c["R_2_F_V_cond0"]++;

	bool L_2_F_cond0 = L_2_F_U_cond0 && L_2_F_V_cond0;
	bool L_2_N_cond0 = L_2_N_U_cond0 && L_2_N_V_cond0;
	bool R_2_N_cond0 = R_2_N_U_cond0 && R_2_N_V_cond0;
	bool R_2_F_cond0 = R_2_F_U_cond0 && R_2_F_V_cond0;

	if (L_2_F_cond0) c["L_2_F_cond0"]++;

	if (L_2_N_cond0) c["L_2_N_cond0"]++;

	if (R_2_N_cond0) c["R_2_N_cond0"]++;

	if (R_2_F_cond0) c["R_2_F_cond0"]++;

	if (L_2_F_cond0 && L_2_N_cond0 && R_2_N_cond0 && R_2_F_cond0)
		c["cond0, U&V, 4RP"]++;

	//--------------------
	
	bool L_2_F_cond1m = (L_2_F_U_cond1m && L_2_F_V_cond0) || (L_2_F_U_cond0 && L_2_F_V_cond1m);
	bool L_2_F_cond1p = (L_2_F_U_cond1p && L_2_F_V_cond0) || (L_2_F_U_cond0 && L_2_F_V_cond1p);
	
	bool L_2_N_cond1m = (L_2_N_U_cond1m && L_2_N_V_cond0) || (L_2_N_U_cond0 && L_2_N_V_cond1m);
	bool L_2_N_cond1p = (L_2_N_U_cond1p && L_2_N_V_cond0) || (L_2_N_U_cond0 && L_2_N_V_cond1p);
	
	bool R_2_N_cond1m = (R_2_N_U_cond1m && R_2_N_V_cond0) || (R_2_N_U_cond0 && R_2_N_V_cond1m);
	bool R_2_N_cond1p = (R_2_N_U_cond1p && R_2_N_V_cond0) || (R_2_N_U_cond0 && R_2_N_V_cond1p);
	
	bool R_2_F_cond1m = (R_2_F_U_cond1m && R_2_F_V_cond0) || (R_2_F_U_cond0 && R_2_F_V_cond1m);
	bool R_2_F_cond1p = (R_2_F_U_cond1p && R_2_F_V_cond0) || (R_2_F_U_cond0 && R_2_F_V_cond1p);

	if (L_2_F_cond1m) c["L_2_F_cond-1"]++;
	if (L_2_F_cond1p) c["L_2_F_cond+1"]++;

	if (L_2_N_cond1m) c["L_2_N_cond-1"]++;
	if (L_2_N_cond1p) c["L_2_N_cond+1"]++;

	if (R_2_N_cond1m) c["R_2_N_cond-1"]++;
	if (R_2_N_cond1p) c["R_2_N_cond+1"]++;

	if (R_2_F_cond1m) c["R_2_F_cond-1"]++;
	if (R_2_F_cond1p) c["R_2_F_cond+1"]++;
	
	//--------------------

	if (L_2_F_recon) c["L_2_F_recon"]++;
	if (L_2_N_recon) c["L_2_N_recon"]++;
	if (R_2_N_recon) c["R_2_N_recon"]++;
	if (R_2_F_recon) c["R_2_F_recon"]++;
	
	bool L_2_F_cond0_AOR = L_2_F_cond0 && L_2_N_recon && R_2_N_recon && R_2_F_recon;
	bool L_2_N_cond0_AOR = L_2_F_recon && L_2_N_cond0 && R_2_N_recon && R_2_F_recon;
	bool R_2_N_cond0_AOR = L_2_F_recon && L_2_N_recon && R_2_N_cond0 && R_2_F_recon;
	bool R_2_F_cond0_AOR = L_2_F_recon && L_2_N_recon && R_2_N_recon && R_2_F_cond0;

	bool L_2_F_cond1m_AOR = L_2_F_cond1m && L_2_N_recon && R_2_N_recon && R_2_F_recon;
	bool L_2_N_cond1m_AOR = L_2_F_recon && L_2_N_cond1m && R_2_N_recon && R_2_F_recon;
	bool R_2_N_cond1m_AOR = L_2_F_recon && L_2_N_recon && R_2_N_cond1m && R_2_F_recon;
	bool R_2_F_cond1m_AOR = L_2_F_recon && L_2_N_recon && R_2_N_recon && R_2_F_cond1m;

	bool L_2_F_cond1p_AOR = L_2_F_cond1p && L_2_N_recon && R_2_N_recon && R_2_F_recon;
	bool L_2_N_cond1p_AOR = L_2_F_recon && L_2_N_cond1p && R_2_N_recon && R_2_F_recon;
	bool R_2_N_cond1p_AOR = L_2_F_recon && L_2_N_recon && R_2_N_cond1p && R_2_F_recon;
	bool R_2_F_cond1p_AOR = L_2_F_recon && L_2_N_recon && R_2_N_recon && R_2_F_cond1p;

	bool rp_any_cond0_AOR = L_2_F_cond0_AOR || L_2_N_cond0_AOR || R_2_N_cond0_AOR || R_2_F_cond0_AOR;
	bool rp_any_cond1m_AOR = L_2_F_cond1m_AOR || L_2_N_cond1m_AOR || R_2_N_cond1m_AOR || R_2_F_cond1m_AOR;
	bool rp_any_cond1p_AOR = L_2_F_cond1p_AOR || L_2_N_cond1p_AOR || R_2_N_cond1p_AOR || R_2_F_cond1p_AOR;

	if (L_2_F_cond0_AOR) c["L_2_F_cond0_AOR"]++;
	if (L_2_N_cond0_AOR) c["L_2_N_cond0_AOR"]++;
	if (R_2_N_cond0_AOR) c["R_2_N_cond0_AOR"]++;
	if (R_2_F_cond0_AOR) c["R_2_F_cond0_AOR"]++;

	if (rp_any_cond0_AOR) c["any_rp_cond0_AOR"]++;
	if (rp_any_cond1m_AOR) c["any_rp_cond-1_AOR"]++;
	if (rp_any_cond1p_AOR) c["any_rp_cond+1_AOR"]++;

	//--------------------

	bool L_2_F_cond0_AOMtr = L_2_F_cond0 && L_2_N_mtr && R_2_N_mtr && R_2_F_mtr;
	bool L_2_N_cond0_AOMtr = L_2_F_mtr && L_2_N_cond0 && R_2_N_mtr && R_2_F_mtr;
	bool R_2_N_cond0_AOMtr = L_2_F_mtr && L_2_N_mtr && R_2_N_cond0 && R_2_F_mtr;
	bool R_2_F_cond0_AOMtr = L_2_F_mtr && L_2_N_mtr && R_2_N_mtr && R_2_F_cond0;

	bool L_2_F_cond1m_AOMtr = L_2_F_cond1m && L_2_N_mtr && R_2_N_mtr && R_2_F_mtr;
	bool L_2_N_cond1m_AOMtr = L_2_F_mtr && L_2_N_cond1m && R_2_N_mtr && R_2_F_mtr;
	bool R_2_N_cond1m_AOMtr = L_2_F_mtr && L_2_N_mtr && R_2_N_cond1m && R_2_F_mtr;
	bool R_2_F_cond1m_AOMtr = L_2_F_mtr && L_2_N_mtr && R_2_N_mtr && R_2_F_cond1m;

	bool L_2_F_cond1p_AOMtr = L_2_F_cond1p && L_2_N_mtr && R_2_N_mtr && R_2_F_mtr;
	bool L_2_N_cond1p_AOMtr = L_2_F_mtr && L_2_N_cond1p && R_2_N_mtr && R_2_F_mtr;
	bool R_2_N_cond1p_AOMtr = L_2_F_mtr && L_2_N_mtr && R_2_N_cond1p && R_2_F_mtr;
	bool R_2_F_cond1p_AOMtr = L_2_F_mtr && L_2_N_mtr && R_2_N_mtr && R_2_F_cond1p;

	bool rp_any_cond0_AOMtr = L_2_F_cond0_AOMtr || L_2_N_cond0_AOMtr || R_2_N_cond0_AOMtr || R_2_F_cond0_AOMtr;
	bool rp_any_cond1m_AOMtr = L_2_F_cond1m_AOMtr || L_2_N_cond1m_AOMtr || R_2_N_cond1m_AOMtr || R_2_F_cond1m_AOMtr;
	bool rp_any_cond1p_AOMtr = L_2_F_cond1p_AOMtr || L_2_N_cond1p_AOMtr || R_2_N_cond1p_AOMtr || R_2_F_cond1p_AOMtr;

	if (L_2_F_cond0_AOMtr) c["L_2_F_cond0_AOMtr"]++;
	if (L_2_N_cond0_AOMtr) c["L_2_N_cond0_AOMtr"]++;
	if (R_2_N_cond0_AOMtr) c["R_2_N_cond0_AOMtr"]++;
	if (R_2_F_cond0_AOMtr) c["R_2_F_cond0_AOMtr"]++;

	if (rp_any_cond0_AOMtr) c["any_rp_cond0_AOMtr"]++;
	if (rp_any_cond1m_AOMtr) c["any_rp_cond-1_AOMtr"]++;
	if (rp_any_cond1p_AOMtr) c["any_rp_cond+1_AOMtr"]++;
}

//----------------------------------------------------------------------------------------------------

void PrintMultiplicityArray(const vector<int> &multiplicities)
{
	for (const auto m : multiplicities)
		printf("%u, ", m);
}

//----------------------------------------------------------------------------------------------------

void TestAlgorithm()
{
	vector<int> multiplicities;
   
	printf("----------------------------\n");
	printf("condition =5, =5, =5, >5, >5; U projection\n");
	printf("----------------------------\n");

	multiplicities = { 0, 5, 0, 8, 0, 5, 0, 5, 0, 9 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pU, {"=5", "=5", "=5", ">5", ">5"}));

	multiplicities = { 0, 1, 0, 1, 0, 0, 0, 1, 0, 1 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pU, {"=5", "=5", "=5", ">5", ">5"}));

	multiplicities = { 0, 5, 0, 5, 0, 5, 0, 5, 0, 5 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pU, {"=5", "=5", "=5", ">5", ">5"}));

	multiplicities = { 0, 9, 0, 9, 0, 5, 0, 5, 0, 5 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pU, {"=5", "=5", "=5", ">5", ">5"}));

	printf("----------------------------\n");
	printf("condition =5, =5, =5, >5, >5; V projection\n");
	printf("----------------------------\n");

	multiplicities = { 0, 5, 0, 8, 0, 5, 0, 5, 0, 9 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pV, {"=5", "=5", "=5", ">5", ">5"}));

	multiplicities = { 5, 0, 8, 0, 5, 0, 5, 0, 9, 0 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pV, {"=5", "=5", "=5", ">5", ">5"}));

	printf("----------------------------\n");
	printf("condition =1, =1, =1, <4, <4; U projection\n");
	printf("----------------------------\n");

	multiplicities = { 0, 5, 0, 8, 0, 5, 0, 5, 0, 9 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pU, {"=1", "=1", "=1", "<4", "<4"}));

	multiplicities = { 0, 3, 0, 1, 0, 1, 0, 1, 0, 2 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pU, {"=1", "=1", "=1", "<4", "<4"}));

	multiplicities = { 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 };
	PrintMultiplicityArray(multiplicities);
	printf("  --> %u\n", FitsMultiplicityCondition(multiplicities, pU, {"=1", "=1", "=1", "<4", "<4"}));
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	/*
	TestAlgorithm();
	return 0;
	*/

	// select input
	vector<string> input_files;

	for (int argi = 1; argi < argc; ++argi)
	{
		if (strcmp(argv[argi], "9009") == 0)
			input_files.push_back("root://eostotem.cern.ch//eos/totem/data/offline/2013/11m/Ntuple/v_4.0/val9009_totem_ntuple.root");

		if (strcmp(argv[argi], "9010") == 0)
			input_files.push_back("root://eostotem.cern.ch//eos/totem/data/offline/2013/11m/Ntuple/v_4.0/val9010_totem_ntuple.root");
	}

	// get input
	TChain *ch = new TChain("TotemNtuple");
	printf(">> input_files\n");
	for (unsigned int i = 0; i < input_files.size(); i++)
	{
		ch->Add(input_files[i].c_str());
		printf("    %s\n", input_files[i].c_str());
	}
	printf("    entries: %llu\n", ch->GetEntries());
	
	// select and link input branches
	ch->SetBranchStatus("*", 0);

	/*
	EventMetaData *metaData = new EventMetaData();
	ch->SetBranchStatus("event_info.*", 1);
	ch->SetBranchAddress("event_info.", &metaData);
	*/

	TriggerData *triggerData = new TriggerData();
	ch->SetBranchStatus("trigger_data.*", 1);
	ch->SetBranchAddress("trigger_data.", &triggerData);

	DiagStruct diag_45b, diag_45t;
	diag_45b.AssignBranches(ch, 25, 21, 120, 124);
	diag_45t.AssignBranches(ch, 24, 20, 121, 125);

	// prepare counters and histograms
	map<string, CounterMap> counters;	// map: diagonal label -> CounterMap

	// loop over events
	for (unsigned int ev = 0; ev < ch->GetEntries(); ev++)
	{
		if ((ev % 100000) == 0)
			printf("%u\n", ev);

		//if (ev >= 100000)
		//	break;

		ch->GetEvent(ev);

		// check bunch number
		if (triggerData->bunch_num != 0)
			continue;
		
		// run analysis
		AnalyzeDiagonal(diag_45b, ev, counters["45b_56t"]);
		AnalyzeDiagonal(diag_45t, ev, counters["45t_56b"]);
	}

	// print results
	for (const auto &p : counters)
	{
		printf(">> %s\n", p.first.c_str());

		for (const auto &e : p.second)
		{
			printf("    %s : %u\n", e.first.c_str(), e.second);
		}
	}

	return 0;
}
