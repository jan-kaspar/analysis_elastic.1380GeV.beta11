import root;
import pad_layout;

string topDir = "../../";

TH2_palette = Gradient(blue, heavygreen, yellow, red, black);

yTicksDef = RightTicks(Step = 50, step = 10);

string dataSets[] = {
	"DS1"
};

string dgns[] = { "45b_56t", "45t_56b" };
string dgn_labs[] = { "45b-56t", "45t-56b" };

//----------------------------------------------------------------------------------------------------

for (int dsi : dataSets.keys)
{
	NewRow();

	for (int dgi : dgns.keys)
	{
		string f = topDir + dataSets[dsi] + "/distributions_" + dgns[dgi] + ".root";

		NewPad("$\ph^*$", "$\th^*\ung{\mu rad}$", axesAbove = false);
		scale(Linear, Linear, Log);
		draw(scale(1., 1e6), RootGetObject(f, "acceptance correction/h_th_vs_phi_after"));
		ylimits(200, 550, Crop);

		/*
		pen ap = cyan;
		for (real y = 0; y <= 140; y += 20)
		{
			if (y >= 80) ap = black;
			xaxis(YEquals(y, false), ap, above=true);
		}
		*/

		AttachLegend(dataSets[dsi]+", "+dgn_labs[dgi], N, N);
	}
}
