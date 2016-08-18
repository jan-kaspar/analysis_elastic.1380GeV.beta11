import root;
import pad_layout;

string topDir = "../../";

//string datasets[] = { "DS1", "DS3" };
string datasets[] = { "DS1" };

string dgns[] = { "45b_56t", "45t_56b" };
string dgn_labs[] = { "45 bot -- 56 top", "45 top -- 56 bot" };

string plot[] = { "p_th_x_diffLR_vs_th_x", "p_th_y_diffLR_vs_th_y", "p_th_y_L_diffNF_vs_th_y_L", "p_th_y_R_diffNF_vs_th_y_R" };
string lab_h[] = { "$\th_x^{*}\ung{\mu rad}$", "$\th_y^{*}\ung{\mu rad}$", "$\th_y^{*L}\ung{\mu rad}$", "$\th_y^{*R}\ung{\mu rad}$" };
string lab_v[] = { "$\De^{R-L}\th_x^{*}\ung{\mu rad}$", "$\De^{R-L}\th_y^{*}\ung{\mu rad}$", "$\De^{F-N}\th_y^{*L}\ung{\mu rad}$", "$\De^{F-N}\th_y^{*R}\ung{\mu rad}$" };
real y_min[] = { -5, -2, -0.1, -0.1};
real y_max[] = { 5, 2, 0.1, 0.1};

for (int dsi : datasets.keys)
{
	for (int dgi : dgns.keys)
	{
		NewRow();

		NewPad(false);	
		label("\vbox{\SetFontSizesXX\hbox{"+datasets[dsi]+"}\hbox{" + dgn_labs[dgi]+"}}");
		
		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgi] + ".root";

		for (int pi : plot.keys)
		{
			NewPad(lab_h[pi], lab_v[pi]);

			string base = "selected - angles/"+plot[pi];

			draw(scale(1e6, 1e6), rGetObj(f, base), "d0,eb", red);

			rObject fit = rGetObj(f, base + "|pol1", error=false);
			if (fit.valid)
			{
				real slope = fit.rExec("GetParameter", 1);
				draw(scale(1e6, 1e6), fit, "", blue+2pt, format("slope = %.4f", slope));
			}

			//ylimits(y_min[pi], y_max[pi], Crop);
			AttachLegend(NW, NW);
		}
	}
}
