//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct DiffPlots
{
	TH1D *h_de_th_x_L, *h_de_th_x_R;
	TH1D *h_de_th_y_L, *h_de_th_y_R;
	TH1D *h_de_th_x, *h_de_th_y, *h_de_th;
	TH1D *h_de_t_x, *h_de_t_y, *h_de_t;
	
	TProfile *p_de_th_x, *p_de_th_y, *p_de_th;
	TProfile *p_de_t_x, *p_de_t_y, *p_de_t;

	DiffPlots();
	void Fill(const Kinematics &act, const Kinematics &nom, double w);
	void Write() const;
};

//----------------------------------------------------------------------------------------------------

DiffPlots::DiffPlots()
{
	unsigned int N_bins = 500;
	h_de_th_x_L = new TH1D("", ";#Delta^{test - true}  #theta_{x}^{L}   (#murad)", N_bins, 0., 0.);
	h_de_th_x_R = new TH1D("", ";#Delta^{test - true}  #theta_{x}^{R}   (#murad)", N_bins, 0., 0.);
	h_de_th_x = new TH1D("", ";#Delta^{test - true}  #theta_{x}   (#murad)", N_bins, 0., 0.);

	h_de_th_y_L = new TH1D("", ";#Delta^{test - true}  #theta_{y}^{L}   (#murad)", N_bins, 0., 0.);
	h_de_th_y_R = new TH1D("", ";#Delta^{test - true}  #theta_{y}^{R}   (#murad)", N_bins, 0., 0.);
	h_de_th_y = new TH1D("", ";#Delta^{test - true}  #theta_{y}   (#murad)", N_bins, 0., 0.);

	h_de_th   = new TH1D(""  , ";#Delta^{test - true}  #theta   (#murad)"    , N_bins, 0., 0.);
	
	h_de_t_x = new TH1D("", ";#Delta^{test - true}  t_{x}   (GeV^{2})", N_bins, 0., 0.);
	h_de_t_y = new TH1D("", ";#Delta^{test - true}  t_{y}   (GeV^{2})", N_bins, 0., 0.);
	h_de_t   = new TH1D(""  , ";#Delta^{test - true}  t   (GeV^{2})"    , N_bins, 0., 0.);

	p_de_th_x = new TProfile("", ";#theta_{x}^{true}   (#murad);mean of #Delta^{test - true} #theta_{x}   (#murad)", N_bins, 0., 0.);
	p_de_th_y = new TProfile("", ";#theta_{y}^{true}   (#murad);mean of #Delta^{test - true} #theta_{y}   (#murad)", N_bins, 0., 0.);
	p_de_th   = new TProfile(""  , ";#theta^{true}   (#murad);mean of #Delta^{test - true} #theta   (#murad)",         N_bins, 0., 0.);

	p_de_t_x = new TProfile("", ";t_{x}^{true}   (GeV^{2});mean of #Delta^{test - true} t_{x}   (GeV^{2})", N_bins, 0., 0.);
	p_de_t_y = new TProfile("", ";t_{y}^{true}   (GeV^{2});mean of #Delta^{test - true} t_{y}   (GeV^{2})", N_bins, 0., 0.);
	p_de_t   = new TProfile(""  , ";t^{true}   (GeV^{2});mean of #Delta^{test - true} t   (GeV^{2})",         N_bins, 0., 0.);
}

//----------------------------------------------------------------------------------------------------

void DiffPlots::Fill(const Kinematics &act, const Kinematics &nom, double w)
{
	h_de_th_x_L->Fill(act.th_x_L - nom.th_x_L, w);
	h_de_th_x_R->Fill(act.th_x_R - nom.th_x_R, w);	
	h_de_th_x->Fill(act.th_x - nom.th_x, w);

	h_de_th_y_L->Fill(act.th_y_L - nom.th_y_L, w);
	h_de_th_y_R->Fill(act.th_y_R - nom.th_y_R, w);	
	h_de_th_y->Fill(act.th_y - nom.th_y, w);

	h_de_th->Fill(act.th - nom.th, w);
	
	h_de_t_x->Fill(act.t_x - nom.t_x, w);
	h_de_t_y->Fill(act.t_y - nom.t_y, w);
	h_de_t->Fill(act.t - nom.t, w);

	p_de_th_x->Fill(nom.th_x, act.th_x - nom.th_x, w);
	p_de_th_y->Fill(nom.th_y, act.th_y - nom.th_y, w);
	p_de_th->Fill(nom.th, act.th - nom.th, w);

	p_de_t_x->Fill(nom.t_x, act.t_x - nom.t_x, w);
	p_de_t_y->Fill(nom.t_y, act.t_y - nom.t_y, w);
	p_de_t->Fill(nom.t, act.t - nom.t, w);
}

//----------------------------------------------------------------------------------------------------

void DiffPlots::Write() const
{
	h_de_th_x_L->Write("h_de_th_x_L");
	h_de_th_x_R->Write("h_de_th_x_R");
	h_de_th_x->Write("h_de_th_x");

	h_de_th_y_L->Write("h_de_th_y_L");
	h_de_th_y_R->Write("h_de_th_y_R");
	h_de_th_y->Write("h_de_th_y");

	h_de_th->Write("h_de_th");
	
	h_de_t_x->Write("h_de_t_x");
	h_de_t_y->Write("h_de_t_y");
	h_de_t->Write("h_de_t");

	p_de_th_x->Write("p_de_th_x");
	p_de_th_y->Write("p_de_th_y");
	p_de_th->Write("p_de_th");

	p_de_t_x->Write("p_de_t_x");
	p_de_t_y->Write("p_de_t_y");
	p_de_t->Write("p_de_t");

	// TODO: add derived plots (de_t/t vs t., ...)
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct CutPlots
{
	map<unsigned int, TH1D *> h_cq;
	map<unsigned int, TH2D *> h2_cq;
	map<unsigned int, TProfile *> p_cq;

	CutPlots(const Analysis &anal);

	void Fill(const CutData &, const Analysis &);
	void Write();
};

//----------------------------------------------------------------------------------------------------

CutPlots::CutPlots(const Analysis &anal)
{
	for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
	//for (unsigned idx = 0; idx < anal.cuts.size(); ++idx)
	{
		//unsigned int ci = anal.cuts[idx];
		char title[100];

		sprintf(title, ";cq%u", ci);
		h_cq[ci] = new TH1D("", title, 150, 0., 0.);

		sprintf(title, ";%s;%s", anal.cqaN[ci].c_str(), anal.cqbN[ci].c_str());
		h2_cq[ci] = new TH2D("", title, 100, 0., 0., 100, 0., 0.);
		p_cq[ci] = new TProfile("", title, 300, 0., 0.);
	}
}

//----------------------------------------------------------------------------------------------------

void CutPlots::Fill(const CutData &cd, const Analysis &anal)
{
	for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
	//for (unsigned idx = 0; idx < anal.cuts.size(); ++idx)
	{
		//unsigned int ci = anal.cuts[idx];
		h_cq[ci]->Fill(cd.cv[ci]);
		h2_cq[ci]->Fill(cd.cqa[ci], cd.cqb[ci]);
		p_cq[ci]->Fill(cd.cqa[ci], cd.cqb[ci]);			
	}
}

//----------------------------------------------------------------------------------------------------

void CutPlots::Write()
{
	TDirectory *baseDir = gDirectory;

	for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
	//for (unsigned idx = 0; idx < anal.cuts.size(); ++idx)
	{
		//unsigned int ci = anal.cuts[idx];

		char buf[100];
		sprintf(buf, "cut %u", ci);
		gDirectory = baseDir->mkdir(buf);

		char name[100];
		sprintf(name, "h_cq%u", ci);	h_cq[ci]->Write(name);
		sprintf(name, "h2_cq%u", ci);	h2_cq[ci]->Write(name);
		sprintf(name, "p_cq%u", ci);	p_cq[ci]->Write(name);
	}

	gDirectory = baseDir;
}
