import root;
import pad_layout;

string topDir = "../../";

//TH2_palette = Gradient(blue, heavygreen, yellow, red, black);
TH2_palette = new pen[] { green, blue, yellow, red, black, magenta };

yTicksDef = RightTicks(Step = 50, step = 10);

xSizeDef = 8cm;

//string dataSets[] = { "DS1",  "DS3" };
string dataSets[] = { "DS1" };

string dgns[] = { "45b_56t", "45t_56b" };
string dgn_labs[] = { "45 bot -- 56 top", "45 top -- 56 top" };

for (int dsi : dataSets.keys)
{
	NewRow();

	for (int dgi : dgns.keys) {
		string f = topDir + dataSets[dsi] + "/distributions_" + dgns[dgi] + ".root";

		NewPad("$\ph^*$", "$\th^*\ung{\mu rad}$", axesAbove = false);
		scale(Linear, Linear, Log);
		draw(scale(1., 1e6), rGetObj(f, "acceptance correction/h_th_vs_phi_after"));

		real x_min = 0, x_max = 0;
		if (dgns[dgi] == "45b_56t") { x_min = 0; x_max = +pi; }
		if (dgns[dgi] == "45t_56b") { x_min = -pi; x_max = +0; }

		real y_min = 150, y_max = 550;

		limits((x_min, y_min), (x_max, y_max), Crop);

		pen ap = black;
		for (real y = y_min+50; y <= y_max; y += 50)
		{
			if (y >= 80)
				ap = black;
			xaxis(YEquals(y, false), ap, above=true);
		}

		AttachLegend(dataSets[dsi]+", "+dgn_labs[dgi], N, N);
	}
}
