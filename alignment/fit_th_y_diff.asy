import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	//"DS1_no_add_alignment",
	"DS1",
	//"DS1_glob_al_45t_56b",
	//"DS1_rel_al_only",
};

string dgns[], dgn_labs[];
dgns.push("45b_56t"); dgn_labs.push("45 bot -- 56 top");
dgns.push("45t_56b"); dgn_labs.push("45 top -- 56 bot");

xSizeDef = 8cm;

xTicksDef = LeftTicks(100., 50.);

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	for (int dgni : dgns.keys)
	{
		string f = "fit_th_y_diff.root";

		real x_min, x_max;
		if (dgns[dgni] == "45b_56t")
		{
			x_min = 0;
			x_max = 500;
			TH1_x_min = +150e-6;
			TH1_x_max = +500e-6;
		} else {
			x_min = -500;
			x_max = 0;
			TH1_x_min = -500e-6;
			TH1_x_max = -150e-6;
		}

		NewRow();

		NewPad(false);
		label("\vbox{\SetFontSizesXX\hbox{" + replace(datasets[dsi], "_", "\_") + "}\hbox{" + dgn_labs[dgni] + "}}");

		string dir = datasets[dsi] + "/" + dgns[dgni] + "/";

		yTicksDef = RightTicks(1., 0.5);

		NewPad("$\th^{*L}_y\ung{\mu rad}$", "$\th_y^{*LF} - \th_y^{*LN} \ung{\mu rad}$");
		draw(scale(1e6, 1e6), RootGetObject(f, dir + "p_th_y_L_diffNF_vs_th_y_L"), "d0,eb", blue);
		RootGetObject(f, dir + "p_th_y_L_diffNF_vs_th_y_L|ff_pol1");
		TF1_x_min = -inf; TF1_x_max = +inf;
		draw(scale(1e6, 1e6), robj, magenta+2pt);
		TF1_x_min = x_min*1e-6; TF1_x_max = x_max*1e-6;
		draw(scale(1e6, 1e6), robj, magenta+dashed);
		limits((x_min, -3), (x_max, +3), Crop);

		NewPad("$\th^{*R}_y\ung{\mu rad}$", "$\th_y^{*RF} - \th_y^{*RN} \ung{\mu rad}$");
		draw(scale(1e6, 1e6), RootGetObject(f, dir + "p_th_y_R_diffNF_vs_th_y_R"), "d0,eb", blue);
		RootGetObject(f, dir + "p_th_y_R_diffNF_vs_th_y_R|ff_pol1");
		TF1_x_min = -inf; TF1_x_max = +inf;
		draw(scale(1e6, 1e6), robj, magenta+2pt);
		TF1_x_min = x_min*1e-6; TF1_x_max = x_max*1e-6;
		draw(scale(1e6, 1e6), robj, magenta+dashed);
		limits((x_min, -3), (x_max, +3), Crop);

		yTicksDef = RightTicks(10., 5.);

		NewPad("$\th^*_y\ung{\mu rad}$", "$\th_y^{*R} - \th_y^{*L} \ung{\mu rad}$");
		draw(scale(1e6, 1e6), RootGetObject(f, dir + "p_th_y_diffLR_vs_th_y"), "d0,eb", blue);
		RootGetObject(f, dir + "p_th_y_diffLR_vs_th_y|ff_pol0");
		TF1_x_min = -inf; TF1_x_max = +inf;
		draw(scale(1e6, 1e6), robj, magenta+2pt);
		TF1_x_min = x_min*1e-6; TF1_x_max = x_max*1e-6;
		draw(scale(1e6, 1e6), robj, magenta+dashed);
		limits((x_min, -40), (x_max, +40), Crop);
	}
}
