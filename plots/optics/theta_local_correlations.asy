import root;
import pad_layout;

string topDir = "../../";

//string datasets[] = { "DS1", "DS3" };
string datasets[] = { "DS1" };

string dgns[] = { "45b_56t", "45t_56b" };
string dgn_labs[] = { "45 bot -- 56 top", "45 top -- 56 bot" };

xSizeDef = 8cm;
ySizeDef = 8cm;

drawGridDef = true;

TH2_z_max = 170;

int y = 0;

for (int dsi : datasets.keys)
{
	y += 2;

	int x = 0;

	NewPad(false, x, y);	
	label("\vbox{\SetFontSizesXX\hbox{"+datasets[dsi]+"}}");

	for (int dgi : dgns.keys)
	{
		x += 1;

		NewPad("$\th^{\rm loc,R}_y\ung{\mu rad}$", "$\th^{\rm loc,L}_y\ung{\mu rad}$", x, y);
		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgi] + ".root";
		draw(scale(1e6, 1e6), rGetObj(f, "optics/p_thl_y_L_vs_thl_y_R"), "eb,d0", StdPen(dgi+1));
		limits((-40, -60), (+40, +60), Crop);
		AttachLegend(dgn_labs[dgi]);
	
		NewPad("$\th^{\rm loc,R}_y\ung{\mu rad}$", "$\th^{\rm loc,L}_y\ung{\mu rad}$", x, y+1);
		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgi] + ".root";
		draw(scale(1e6, 1e6), rGetObj(f, "optics/h_thl_y_L_vs_thl_y_R"), "def");
		limits((-40, -60), (+40, +60), Crop);
		AttachLegend(dgn_labs[dgi]);
	}

	x += 1;

	NewPad("$y^{LF}\ung{mm}$", "$\th_y^{\rm loc,L}\ung{\mu rad}$", x, y);
	for (int dgi : dgns.keys)
	{
		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgi] + ".root";
		draw(scale(1e3, 1e6), rGetObj(f, "optics/p_thl_y_L_vs_y_LF"), "eb,d0", StdPen(dgi+1), dgn_labs[dgi]);
	}
	limits((-10, -50), (+10, +50), Crop);

	NewPad("$y^{LF}\ung{mm}$", "$\th_y^{\rm loc,L}\ung{\mu rad}$", x, y+1);
	for (int dgi : dgns.keys)
	{
		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgi] + ".root";
		draw(scale(1e3, 1e6), rGetObj(f, "optics/h_thl_y_L_vs_y_LF"), (dgi == 0) ? "p,d0,bar" : "p,d0");
	}
	limits((-10, -50), (+10, +50), Crop);
	
	x += 1;

	NewPad("$y^{RF}\ung{mm}$", "$\th_y^{\rm loc,R}\ung{\mu rad}$", x, y);
	for (int dgi : dgns.keys)
	{
		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgi] + ".root";
		draw(scale(1e3, 1e6), rGetObj(f, "optics/p_thl_y_R_vs_y_RF"), "eb,d0", StdPen(dgi+1), dgn_labs[dgi]);
	}
	limits((-10, -50), (+10, +50), Crop);

	NewPad("$y^{RF}\ung{mm}$", "$\th_y^{\rm loc,R}\ung{\mu rad}$", x, y+1);
	for (int dgi : dgns.keys)
	{
		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgi] + ".root";
		draw(scale(1e3, 1e6), rGetObj(f, "optics/h_thl_y_R_vs_y_RF"), (dgi == 0) ? "p,d0,bar" : "p,d0");
	}
	limits((-10, -50), (+10, +50), Crop);

}
