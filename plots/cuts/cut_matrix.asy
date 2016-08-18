import root;
import pad_layout;

string topDir = "../../";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

string datasets[] = { "DS1" };

// TODO
//string dgns[] = { "45b_56t", "45t_56b" };
//string dgns[] = { "45b_56t" };
string dgns[] = { "45t_56b" };

int cuts[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

real scale_x[] = { 1e6, 1e6, 1e6, 1e6, 1e0, 1e0, 1e6, 1e6 };
real scale_y[] = { 1e6, 1e6, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0 };

string label_x[] = { "$\th_x^{*R}\ung{\mu rad}$", "$\th_y^{*R}\ung{\mu rad}$", "$\th_x^{*R}\ung{\mu rad}$", "$\th_x^{*L}\ung{\mu rad}$", "$y^{R,N}\ung{mm}$", "$y^{L,N}\ung{mm}$", "$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$" };
string label_y[] = { "$\th_x^{*L}\ung{\mu rad}$", "$\th_y^{*L}\ung{\mu rad}$", "$x^{*R}\ung{mm}$", "$x^{*L}\ung{mm}$", "$y^{R,F} - y^{R,N}\ung{mm}$", "$y^{L,F} - y^{L,N}\ung{mm}$", "$\De^{R-L} x^*\ung{mm}$", "$\De^{R-L} y^*\ung{mm}$" };
string label_cut[] = { "$\De^{R-L} \th_x^{*}\ung{\mu rad}$", "$\De^{R-L} \th_y^{*}\ung{\mu rad}$", "$x^{*R}\ung{mm}$", "$x^{*L}\ung{mm}$", "$cq5$", "$cq6$", "$cq7$", "$cq8$" };

real lim_x_low[] = { -1000, -1000, -1000, -1000, -15, -15, -600, -600 };
real lim_x_high[] = { +1000, +1000, +1000, +1000, +15, +15, +600, +600 };

real lim_y_low[] = { -1000, -1000, -0.8, -0.8, -0.5, -0.5, -1.5, -4 };
real lim_y_high[] = { +1000, +1000, +0.8, +0.8, +0.5, +0.5, +1.5, +4 };

string cut_combinations[] = { "no_cuts", "1", "2", "3", "4", "5", "6", "7", "1,2", "1,2,7", "1,2,5,6,7", "1,2,3,4,5,6,7" };

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	string dataset = datasets[dsi];

	write(dataset);

	for (int dgi : dgns.keys)
	{
		string dgn = dgns[dgi];
	
		write("    " + dgn);
		
		NewPad(false);
		for (int ci : cuts.keys)
		{
			NewPad(false);
			label("{\SetFontSizesXX " + format("cut %i", cuts[ci]) + "}");
		}

		for (int cci : cut_combinations.keys)
		{
			write("        " + cut_combinations[cci]);

			NewRow();

			NewPad(false);
			label("\vbox{\SetFontSizesXX\hbox{" + replace(cut_combinations[cci], "_", " ") +"}}");

			string dir = (cut_combinations[cci] == "no_cuts") ? "no_cuts" : "cuts:" + cut_combinations[cci];

			string f = topDir + dataset+"/background_study/" + dir + "/distributions_" + dgn + ".root";
			
			
			for (int ci : cuts.keys)
			{
				int cut = cuts[ci];
				int idx = cut - 1;

				NewPad(label_x[idx], label_y[idx]);
				scale(Linear, Linear, Log);
				string objC = format("elastic cuts/cut %i", cut) + format("/plot_after_cq%i", cut);
				draw(scale(scale_x[idx], scale_y[idx]), rGetObj(f, objC+"#0"), "p,d0,bar");
				draw(scale(scale_x[idx], scale_y[idx]), rGetObj(f, objC+"#1"));
				draw(scale(scale_x[idx], scale_y[idx]), rGetObj(f, objC+"#2"));
				limits((-lim_x_high[idx], -lim_y_high[idx]), (-lim_x_low[idx], -lim_y_low[idx]), Crop);
				
				AttachLegend();
			}
		}
	
		GShipout("cut_matrix_" + dataset + "_" + dgn);
	}
}
