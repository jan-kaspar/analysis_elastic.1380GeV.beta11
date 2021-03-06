import root;
import pad_layout;

string topDir = "./";

//string dataset = "DS1_rel_al_only";

string diagonal = "45b_56t";

string rp = "R_N";

string steps[];
steps.push("sh = -20.00");
steps.push("sh = -10.00");
steps.push("sh = -0.00");
steps.push("sh = 5.00");
steps.push("sh = 10.00");
steps.push("sh = 20.00");

real x_min = 160;
real x_max = 480;

xSizeDef = 8cm;

xTicksDef = LeftTicks(50., 10.);

//----------------------------------------------------------------------------------------------------

NewPad(false);
//label(replace("\vbox{\hbox{" + dataset + "}\hbox{" + diagonal + "}\hbox{" + rp + "}}", "_", "\_"));
label(replace("\vbox{\hbox{" + diagonal + "}\hbox{" + rp + "}}", "_", "\_"));

//----------------------------------------------------------------------------------------------------

for (int sti : steps.keys)
{
	NewRow();

	NewPad(false);
	label(steps[sti]);

	string f = topDir + "match_th_y_scale_fcn.root";
	string base = diagonal + "/" + rp + "/" + steps[sti];

	NewPad("$\th^*_x$ (blue), $\th^*_y$ (red) $\ung{\mu rad}$", "entries");
	draw(scale(1e6, 1), RootGetObject(f, base + "/h_cmp|g_h_th_x"), "l", blue);
	draw(scale(1e6, 1), RootGetObject(f, base + "/h_cmp|h_th_y"), "vl", red);
	xlimits(x_min, x_max, Crop);

	NewPad("$\th^*_x$ (blue), $\th^*_y$ (red) $\ung{\mu rad}$", "cumulative distribution");
	draw(scale(1e6, 1), RootGetObject(f, base + "/ch_cmp|g_ch_th_x"), "l", blue);
	draw(scale(1e6, 1), RootGetObject(f, base + "/ch_cmp|g_ch_th_y"), "l", red);
	xlimits(x_min, x_max, Crop);

	NewPad("$\th^*_{x,y} \ung{\mu rad}$", "cumulative distribution difference");
	draw(scale(1e6, 1), RootGetObject(f, base + "/g_kol_diff_graph"), "l", heavygreen);
	xlimits(x_min, x_max, Crop);
}

GShipout(vSkip=1mm);
