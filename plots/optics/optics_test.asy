import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS1", "DS3" };

string dgns[] = { "45b_56t", "45t_56b" };
string dgn_labs[] = { "far dgn", "close dgn" };

string proj[] = { "1", "2" };
string lab_h[] = { "$\th_x^{*R}\ung{\mu rad}$", "$\th_y^{*R}\ung{\mu rad}$" };
string lab_v[] = { "$\th_x^{*L}\ung{\mu rad}$", "$\th_y^{*L}\ung{\mu rad}$" };

for (int dsi : datasets.keys) {
	for (int dgi : dgns.keys) {
		NewRow();

		NewPad(false);	
		label(datasets[dsi]+", " + dgn_labs[dgi]);
		
		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgi] + ".root";

		for (int pi : proj.keys) {
			NewPad(lab_h[pi], lab_v[pi]);

			draw(scale(1e6, 1e6), rGetObj(f, "elastic cuts/cut "+proj[pi]+"/p_cq"+proj[pi]), "d0,eb", red);

			rObject fit = rGetObj(f, "elastic cuts/cut "+proj[pi]+"/p_cq"+proj[pi]+"|pol1");
			real slope = fit.rExec("GetParameter", 1);
			draw(scale(1e6, 1e6), fit, "", blue, format("slope = %.3f", slope));

			AttachLegend(NW, NW);
		}
	}
}
