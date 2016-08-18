#include "../common_definitions.h"
#include "parameters.h"
#include "../common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile.h"

#include <vector>
#include <map>
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	Init(argv[1]);
	if (diagonal != dCombined)
		return rcIncompatibleDiagonal;

	vector<string> diagonals;
	diagonals.push_back("45b_56t");
	diagonals.push_back("45t_56b");

	vector<TFile *> files;
	for (unsigned int i = 0; i < diagonals.size(); i++)
		files.push_back(new TFile((string("distributions_" + diagonals[i] + ".root").c_str())));

	map<string, string> list;
	list["th_x"] = "acceptance correction/h_th_x_after";
	list["th_y"] = "acceptance correction/h_th_y_after";
	list["th_y_before"] = "acceptance correction/h_th_y_before";
	list["t"] = "acceptance correction/t, ub/h_t_ub_after";
	list["timestamp_dgn"] = "metadata/h_timestamp_dgn";
	list["timestamp_sel"] = "metadata/h_timestamp_sel";

	TFile *outF = new TFile("compare_diagonals.root", "recreate");

	int colors[2] = {2, 4};

	for (map<string, string>::iterator it = list.begin(); it != list.end(); ++it) {
		TCanvas *c = new TCanvas(it->first.c_str());
		c->SetLogy(1);
		vector<TH1D *> hists;
		for (unsigned int d = 0; d < diagonals.size(); d++) {
			string obj = it->second;
			if (it->first.compare("th_y") == 0 && d == 1)
				obj += "_flipped";
			if (it->first.compare("th_y_before") == 0 && d == 1)
				obj += "_flipped";
			TH1D *h = (TH1D *) files[d]->Get(obj.c_str());
			if (!h) {
				printf("ERROR: can't load object `%s'.\n", obj.c_str());
				continue;
			}

			h->Sumw2();
			hists.push_back(h);
			h->SetName(diagonals[d].c_str());
			h->SetLineColor(colors[d]);
			h->Draw((d == 0) ? "" : "sames");
		}
		c->Write();

		TH1D *h_ratio = new TH1D(*hists[0]);
		h_ratio->SetName((it->first + ":ratio").c_str());
		h_ratio->Divide(hists[1]);
		h_ratio->SetLineColor(1);
		h_ratio->Write();
	}

	//--------------------------------------------------

	// th_y_L vs. th_y_R from both diagonals
	TH2D *h1 = (TH2D *) files[0]->Get("selected - angles/h_th_y_L_vs_th_y_R");
	TH2D *h2 = (TH2D *) files[1]->Get("selected - angles/h_th_y_L_vs_th_y_R");

	if (!h1 || !h2) {
		printf("ERROR: can't load histograms `%s'\n", "h_th_y_L_vs_th_y_R");
	} else {
		TH2D *h12 = new TH2D(*h1);
		for (int ix = 1; ix <= h12->GetNbinsX(); ix++) {
			for (int iy = 1; iy <= h12->GetNbinsX(); iy++) {
				double v1 = h1->GetBinContent(ix, iy);
				double v2 = h2->GetBinContent(ix, iy);
				h12->SetBinContent(ix, iy, max(v1, v2));
			}
		}
		h12->Write();
	}
	
	//--------------------------------------------------
	
	TGraph *g1 = (TGraph *) files[0]->Get("selected - angles/g_th_y_L_vs_th_y_R");
	TGraph *g2 = (TGraph *) files[1]->Get("selected - angles/g_th_y_L_vs_th_y_R");

	if (!g1 || !g2) {
		printf("ERROR: can't load graphs `%s'\n", "g_th_y_L_vs_th_y_R");
	} else {
		TGraph *g12cut = new TGraph(*g1);
		g12cut->Set(0);
		
		for (int i = 0; i < g1->GetN(); i++) {
			double x, y;
			g1->GetPoint(i, x, y);
			double th = (x+y)/2.;
			//if (fabs(th) > 30E-6 && fabs(th) < 105E-6)
				g12cut->SetPoint(g12cut->GetN(), x, y);
		}
	
		for (int i = 0; i < g2->GetN(); i++) {
			double x, y;
			g2->GetPoint(i, x, y);
			double th = (x+y)/2.;
			//if (fabs(th) > 30E-6 && fabs(th) < 105E-6)
				g12cut->SetPoint(g12cut->GetN(), x, y);
		}
		g12cut->Fit("pol1");
		g12cut->Write();
	}

	//--------------------------------------------------

	// combine profiles
	map<string, string> profile_list;
	
	profile_list["gp_x_vs_y_L_F"] = "selected - hits/p_x_vs_y_L_F";
	profile_list["gp_x_vs_y_L_N"] = "selected - hits/p_x_vs_y_L_N";
	profile_list["gp_x_vs_y_R_N"] = "selected - hits/p_x_vs_y_R_N";
	profile_list["gp_x_vs_y_R_F"] = "selected - hits/p_x_vs_y_R_F";
	
	profile_list["gp_x_vs_y_L_F_noal"] = "selected - hits/p_x_vs_y_L_F_noal";
	profile_list["gp_x_vs_y_L_N_noal"] = "selected - hits/p_x_vs_y_L_N_noal";
	profile_list["gp_x_vs_y_R_N_noal"] = "selected - hits/p_x_vs_y_R_N_noal";
	profile_list["gp_x_vs_y_R_F_noal"] = "selected - hits/p_x_vs_y_R_F_noal";

	profile_list["gp_th_x_L_vs_th_y_L"] = "selected - angles/p_th_x_L_vs_th_y_L";
	profile_list["gp_th_x_R_vs_th_y_R"] = "selected - angles/p_th_x_R_vs_th_y_R";

	for (map<string, string>::iterator pit = profile_list.begin(); pit != profile_list.end(); ++pit) {
		TGraphErrors *gp = new TGraphErrors(); gp->SetName(pit->first.c_str());

		TProfile *p_0 = (TProfile *) files[0]->Get(pit->second.c_str());
		TProfile *p_1 = (TProfile *) files[1]->Get(pit->second.c_str());

		if (!p_0 || !p_1) {
			printf("ERROR: can't load profiles `%s'\n", pit->second.c_str());
			continue;
		}

		char buf[100];
		sprintf(buf, ";%s;%s", p_0->GetXaxis()->GetTitle(), p_0->GetYaxis()->GetTitle());
		gp->SetTitle(buf);

		for (int bi = 1; bi <= p_0->GetNbinsX(); bi++)
			if (p_0->GetBinError(bi) > 0) {
				int idx = gp->GetN();
				gp->SetPoint(idx, p_0->GetBinCenter(bi), p_0->GetBinContent(bi));
				gp->SetPointError(idx, p_0->GetBinWidth(bi)/2, p_0->GetBinError(bi));
			}
		
		for (int bi = 1; bi <= p_1->GetNbinsX(); bi++)
			if (p_1->GetBinError(bi) > 0) {
				int idx = gp->GetN();
				gp->SetPoint(idx, p_1->GetBinCenter(bi), p_1->GetBinContent(bi));
				gp->SetPointError(idx, p_1->GetBinWidth(bi)/2, p_1->GetBinError(bi));
			}
	
		//gp->Fit("pol1");
		gp->Write();
	}


	delete outF;
	return 0;
}
