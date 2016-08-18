#include "TF1.h"
#include "TRandom2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

namespace CF {

double t_middle = 0.;

unsigned long N_anal_cf = 0;

void MakeFit(TH1D *input, TH1D *corr, TF1 *ff, TH1D *guide)
{
	// correct the input with current correction
	TH1D *input_corr = new TH1D(*input);
	input_corr->Multiply(corr);
	
	input_corr->Fit(ff, "Q",  "", anal.t_min, anal.t_max);
	input_corr->Fit(ff, "Q",  "", anal.t_min, anal.t_max);
	input_corr->Fit(ff, "IQ",  "", anal.t_min, anal.t_max);
	input_corr->Fit(ff, "I",  "", anal.t_min, anal.t_max);

	printf("residual chi^2 = %.1E, number of points = %i, chi^2/npx = %.2f\n",
		ff->GetChisquare(), ff->GetNpx(), ff->GetChisquare() / ff->GetNpx());

	input_corr->SetName("input_corr");
	input_corr->Write();
	
	ff->Write();
}

//----------------------------------------------------------------------------------------------------

void CalculateCorrection(TF1 *ff, TH1D *corr)
{
	TH1D *h_t, *h_t_p;

	Profile *p_t, *p_t_p;
	if (corr->GetXaxis()->GetXbins()->fN > 0) {
		p_t = new Profile("h_t", ";t", corr->GetNbinsX(), corr->GetXaxis()->GetXbins()->fArray);
		p_t_p = new Profile("h_t_p", ";t'", corr->GetNbinsX(), corr->GetXaxis()->GetXbins()->fArray);
	} else {
		p_t = new Profile("h_t", ";t", corr->GetNbinsX(), corr->GetXaxis()->GetXmin(), corr->GetXaxis()->GetXmax());
		p_t_p = new Profile("h_t_p", ";t'", corr->GetNbinsX(), corr->GetXaxis()->GetXmin(), corr->GetXaxis()->GetXmax());
	}
	
	// make simulation
	for (unsigned long ev = 0; ev < N_anal_cf; ++ev) {
		// TODO: update using new classes - Kinematics etc.

		double t = gRandom->Rndm() * (anal.t_max_full - anal.t_min_full) + anal.t_min_full;
		double phi = gRandom->Rndm() * M_PI; // just upper half
		double w = ff->Eval(t) * (anal.t_max_full - anal.t_min_full);

		double th = sqrt(t) / env.p;
		double th_x = th * cos(phi);
		double th_y = th * sin(phi);

		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// TODO: since single and double-arm theta_x resolutions can be radically different, simulate just the double-arm one (the only needed!)
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		double th_x_L_sm = th_x + gRandom->Gaus() * anal.si_th_x_1arm_L;
		double th_x_R_sm = th_x + gRandom->Gaus() * anal.si_th_x_1arm_R;
		double th_y_L_sm = th_y + gRandom->Gaus() * anal.si_th_y_1arm;
		double th_y_R_sm = th_y + gRandom->Gaus() * anal.si_th_y_1arm;

		double th_x_sm = (th_x_L_sm + th_x_R_sm) / 2.;
		double th_y_sm = (th_y_L_sm + th_y_R_sm) / 2.;

		double th_sm = sqrt(th_x_sm*th_x_sm + th_y_sm*th_y_sm);
		double t_sm = env.p*env.p * (th_x_sm*th_x_sm + th_y_sm*th_y_sm);
	
		p_t->Fill(t, w);

		double phi_corr = 0., div_corr = 0.;
		bool skip = true;
		// TODO:
		/*
		bool skip = CalculateAcceptanceCorrections(sign,
			th_x_L_sm, th_x_R_sm, th_x_sm, si_th_x_os,
			th_y_L_sm, th_y_R_sm, th_y_sm, si_th_y_os,
			th_sm, phi_corr, div_corr);
		*/

		if (!skip)
			p_t_p->Fill(t_sm, w * div_corr * phi_corr / 2.); // only upper half in generation
	}
	
	h_t = p_t->GetHist(N_anal_cf);
	h_t_p = p_t_p->GetHist(N_anal_cf);

	h_t->SetLineColor(4); h_t->Write();
	h_t_p->SetLineColor(2); h_t_p->Write();

	// extract correction histogram
	for (int i = 1; i <= h_t->GetNbinsX(); i++) {
		double t = h_t->GetBinCenter(i);
		if (t < anal.t_min || t > anal.t_max) {
			corr->SetBinContent(i, 0.);
			corr->SetBinError(i, 0.);
			continue;
		}

		double c = h_t->GetBinContent(i), e = h_t->GetBinError(i);
		double c_p = h_t_p->GetBinContent(i), e_p = h_t_p->GetBinError(i);

		double v = (c_p > 0.) ? c / c_p : 0.;
		double u = (c_p > 0.) ? sqrt(e*e/c_p/c_p + c*c/c_p/c_p/c_p/c_p * e_p*e_p) : 0.;

		corr->SetBinContent(i, v);
		corr->SetBinError(i, u);
	}

	corr->Write();
}

//----------------------------------------------------------------------------------------------------

vector<TF1 *> parameterizations;

TF1* InitFitFunctions(double norm)
{
	printf(">> CF::InitFitFunctions(%E)\n", norm);

	TF1 *ff;
	
	/*
	ff = new TF1("exp2+exp2", "[0]*exp([1]*x + [2]*x*x) + [3]*exp([4]*x + [5]*x*x)");
	ff->SetParameters(norm*1E2, -20., 0., norm/1E1, -3., 0.);
	ff->SetRange(anal.t_min_full, anal.t_max_full);
	parameterizations.push_back(ff);
	*/
	
	ff = new TF1("exp3+exp2", "[0]*exp([1]*x + [2]*x*x + [3]*x*x*x) + [4]*exp([5]*x + [6]*x*x)");
	ff->SetParameters(norm*1E2, -20., 0., 0., norm/1E1, -3., 0.);
	ff->SetRange(anal.t_min_full, anal.t_max_full);
	parameterizations.push_back(ff);

	/*
	ff = new TF1("exp2+exp4", "[0]*exp([1]*x + [2]*x*x) + [3]*exp([4]*x + [5]*x*x + [6]*x*x*x + [7]*x*x*x*x)");
	ff->SetParameters(norm*1E2, -20., 0., norm/1E1, -3., 0., 0., 0.);
	ff->SetRange(anal.t_min_full, anal.t_max_full);
	parameterizations.push_back(ff);
	*/
	
	ff = new TF1("exp3+exp4", "[0]*exp([1]*x + [2]*x*x + [3]*x*x*x) + [4]*exp([5]*x + [6]*x*x + [7]*x*x*x + [8]*x*x*x*x)");
	ff->SetParameters(norm*1E2, -20., 0., 0., norm/1E1, -3., 0., 0., 0.);
	ff->SetRange(anal.t_min_full, anal.t_max_full);
	parameterizations.push_back(ff);

	return ff;
}

//----------------------------------------------------------------------------------------------------

void DoUnfolding(TH1D *input, TH1D *guide, TF1 *parameterization)
{
	printf(">> CF::DoUnfolding\n");

	// make initial correction histogram
	TH1D *corr = new TH1D(*input);
	corr->SetName("corr");
	for (int i = 1; i <= corr->GetNbinsX(); i++) {
		corr->SetBinContent(i, 1.);
		corr->SetBinError(i, 0.);
	}

	// initialize fit function
	TF1 *ff = new TF1(*parameterization);
	ff->SetName("ff");

	TDirectory *baseDir = gDirectory;

	// run iterations
	for (unsigned int it = 1; it <= 3; it++) {
		printf("* iteration: %i\n", it);
	
		char buf[100];
		sprintf(buf, "iteration %i", it);
		gDirectory = baseDir->mkdir(buf);

		MakeFit(input, corr, ff, guide);
		CalculateCorrection(ff, corr);
	}

	gDirectory = baseDir;

	// save final correction
	corr->SetName("corr_final");
	corr->Write();

	// save unsmeared result
	TH1D *output = new TH1D(*input);
	output->SetName("output");
	output->SetLineColor(8);
	output->Multiply(corr);
	output->Write();

	// distribution comparison
	TCanvas *c = new TCanvas("dist cmp");
	c->SetLogy(1);
	input->SetLineColor(1);
	input->Draw();
	output->Draw("same");
	c->Write();
}

} // namespace
