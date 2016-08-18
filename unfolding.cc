using namespace std;

#include "../common_definitions.h"
#include "parameters.h"
#include "../common.h"
#include "../Profile.h"

#include "../unfolding_CF.h"
//#include "../unfolding_GR.h"

#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom2.h"


//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal settings
	Init(argv[1]);
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	// get input data
	string inFn = string("distributions_") + argv[1] + ".root";	
	TFile *inF = new TFile(inFn.c_str());
	
	TH1D *h_t_ub = (TH1D *) inF->Get("acceptance correction/h_t_ub_after");
	TH1D *h_t_eb = (TH1D *) inF->Get("acceptance correction/h_t_eb_after");
	
	// get EPL 95 data as hight-|t| guide
	/*
	TFile *oldF = new TFile("../publication1.root");
	TH1D *h_t_old = (TH1D *) oldF->Get("h_dsdt");
	h_t_old->Sumw2();
	double sf = full_norm_corr / L_int_eff;
	h_t_old->Scale(1. / sf);
	*/

	// prepare output
	string outFn = string("unfolding_") + argv[1] + ".root";	

	/*
	string outFn = string("unfolding_x+1,y-1") + argv[1] + ".root";	
	si_th_x_os += si_th_x_os_unc;
	si_th_y_os -= si_th_y_os_unc;
	*/

	/*
	string outFn = string("unfolding_x+1,y+1") + argv[1] + ".root";	
	si_th_x_os += si_th_x_os_unc;
	si_th_y_os += si_th_y_os_unc;
	*/

	TFile *outF = new TFile(outFn.c_str(), "recreate");

	printf("si_th_x_1arm_L = %E\n", anal.si_th_x_1arm_L);
	printf("si_th_x_1arm_R = %E\n", anal.si_th_x_1arm_R);
	printf("si_th_y_1arm = %E\n", anal.si_th_y_1arm);

	// set ranges
	//CF::t_middle = 1.86;

	// run CF unfolding
	CF::N_anal_cf = 1E8;
	CF::InitFitFunctions(h_t_eb->GetEntries());

	TDirectory *ubDir = outF->mkdir("cf,ub");
	TDirectory *ebDir = outF->mkdir("cf,eb");

	for (unsigned int pi = 0; pi < CF::parameterizations.size(); ++pi) {
		//gDirectory = ubDir->mkdir(CF::parameterizations[pi]->GetName());
		//CF::DoUnfolding(h_t_ub, NULL, CF::parameterizations[pi]);
	
		gDirectory = ebDir->mkdir(CF::parameterizations[pi]->GetName());
		CF::DoUnfolding(h_t_eb, NULL, CF::parameterizations[pi]);
	}
	
	/*
	// get reference unfolding matrix
	TFile *rF = new TFile("mc_ref.root");
	const TMatrixD &ub_M = GR::GetMatrix((TH2D *) rF->Get("uni/h_res_mat2"));
	const TMatrixD &eb_M = GR::GetMatrix((TH2D *) rF->Get("exp/h_res_mat2"));

	// run GR unfolding
	alpha_gr = 1E-2;
	gDirectory = outF->mkdir("gr,ub");
	GR::DoUnfolding(h_t_ub, ub_M);
	
	alpha_gr = 1E-3;
	gDirectory = outF->mkdir("gr,eb");
	GR::DoUnfolding(h_t_eb, eb_M);
	*/
}
