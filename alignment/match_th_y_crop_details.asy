import root;
import pad_layout;

string topDir = "./";

//string dataset = "DS1_rel_al_only/";
string dataset = "";

string diagonal = "45b_56t";

string rp = "L_F";

string dirs[];
//dirs.push("sh = -50.00");
dirs.push("sh = -0.00");
dirs.push("sh = 10.00");
dirs.push("sh = 20.00");
dirs.push("sh = 30.00");
dirs.push("sh = 40.00");
dirs.push("sh = 50.00");

//----------------------------------------------------------------------------------------------------

for (int diri : dirs.keys)
{
	NewPad(false);
	label(dirs[diri]);
}

NewRow();

for (int diri : dirs.keys)
{
	NewPad("$\th^*_x$ (blue), $\th^*_y$ (red) $\ung{\mu rad}$", "entries");

	string f = topDir + "match_th_y_crop.root";
	string base = dataset + diagonal + "/" + rp + "/" + dirs[diri];
	draw(scale(1e6, 1), RootGetObject(f, base + "/h_cmp|h_th_x"), "vl", blue);
	draw(scale(1e6, 1), RootGetObject(f, base + "/h_cmp|h_th_y"), "vl", red);

	xlimits(150, 450, Crop);
}

NewRow();

for (int diri : dirs.keys)
{
	NewPad("$\th^*_x$ (blue), $\th^*_y$ (red) $\ung{\mu rad}$", "cumulative distribution");

	string f = topDir + "match_th_y_crop.root";
	string base = dataset + diagonal + "/" + rp + "/" + dirs[diri];
	draw(scale(1e6, 1), RootGetObject(f, base + "/ch_cmp|g_ch_th_x"), "l", blue);
	draw(scale(1e6, 1), RootGetObject(f, base + "/ch_cmp|g_ch_th_y"), "l", red);

	xlimits(150, 450, Crop);
}

GShipout(vSkip=1mm);
