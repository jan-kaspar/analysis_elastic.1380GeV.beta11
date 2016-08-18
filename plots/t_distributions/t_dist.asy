import root;
import pad_layout;

string datasets[] = { "DS1", "DS3" };
//string datasets[] = { "DS1" };

string topDir = "../../";

real shifts[] = { 0, -2 };

string dgns[] = { "45b_56t", "45t_56b" };

string binning = "eb";

int idx = 0;

NewPad("$|t|\ung{GeV^2}$", "$\d\si_{\rm el} / \d t \ung{arbitrary\ units}$", 12cm, xTicks=LeftTicks(0.1, 0.05));
scale(Linear, Log);

for (int dsi : datasets.keys) {
	string dataset = datasets[dsi];
	
	for (int dgi : dgns.keys) {
		string dgn = dgns[dgi];

		string f = topDir + dataset+"/distributions_"+dgn+".root";

		pen p = StdPen(idx);
		++idx;

		draw(shift(0, shifts[dsi]), rGetObj(f, "normalization/"+binning+"/h_t_normalized"), "eb", p, replace(dataset + ", " + dgn, "_", "--"));
	}
}


draw(Label("dip expected", 0, N), (0.62, -1)--(0.62, -4), EndArrow);

limits((0, 1e-5), (0.7, 1e3), Crop);
AttachLegend();
