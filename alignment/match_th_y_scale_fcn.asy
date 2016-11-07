import root;
import pad_layout;

string topDir = "./";

//string dataset = "DS1_rel_al_only";

string diagonals[], dgn_labels[];
diagonals.push("45b_56t"); dgn_labels.push("45 bot -- 56 top");
diagonals.push("45t_56b"); dgn_labels.push("45 top -- 56 bot");

string rps[], rp_labels[];
rps.push("L_F"); rp_labels.push("L-F");
rps.push("L_N"); rp_labels.push("L-N");
rps.push("R_N"); rp_labels.push("R-N");
rps.push("R_F"); rp_labels.push("R-F");

xSizeDef = 8cm;

xTicksDef = LeftTicks(10., 2.);
//xTicksDef = LeftTicks(20., 10.);

//----------------------------------------------------------------------------------------------------

for (int dgni : diagonals.keys)
{
	NewPad("$\th^*_y$ shift $\ung{\mu rad}$", "test metric");
	
	for (int rpi : rps.keys)
	{
		string f = topDir + "match_th_y_scale_fcn.root";
		string base = diagonals[dgni] + "/" + rps[rpi];

		pen p = StdPen(rpi);

		draw(scale(1e6, 1), RootGetObject(f, base + "/g_kol_meas"), "l", p, rp_labels[rpi]);
	}

	xlimits(-40, +40, Crop);

	AttachLegend(dgn_labels[dgni]);
}
