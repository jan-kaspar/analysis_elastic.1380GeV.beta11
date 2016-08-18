import root;
import pad_layout;

string topDir = "../../";

//string datasets[] = { "DS1", "DS3" };
string datasets[] = { "DS1" };

string arms[] = { "L", "R" };

xSizeDef = 8cm;
//ySizeDef = 8cm;

xTicksDef = LeftTicks(5., 1.);

//----------------------------------------------------------------------------------------------------

void PlotGroup(string confs[], real y_max=0)
{
	for (int dsi : datasets.keys)
	{
		for (int ai : arms.keys)
		{
			NewPad("multiplicity threshold", "event fraction\ung{\%}");
			
			int idx = 0;
			for (int ci : confs.keys)
			{
				int scp = find(confs[ci], ";");
				string c = substr(confs[ci], 0, scp);
				string v = substr(confs[ci], scp+1);

				string f = topDir+datasets[dsi]+"/showers_combined.root";

				// get normalisation
				rObject h = rGetObj(f, "L/N_T/h_avg_mult");
				real N = h.rExec("GetEntries");

				string p = arms[ai]+"/"+c+"/h_" + v;
				string l = c + ", " + v;
				draw(scale(1., 100./N), rGetObj(f, p), "vl", StdPen(idx), replace(l, "_", "\_"));
				++idx;
			}

			if (y_max > 0)
				ylimits(0, y_max, Crop);

			AttachLegend();
		}
	}
}

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	for (int ai : arms.keys)
	{
		NewPad(false);
		label("{\SetFontSizesXX " + datasets[dsi] + ", arm " + arms[ai] + "}");
	}
}

//----------------------------------------------------------------------------------------------------
NewRow();

string confs[] = {
	"N_T,F_T;inc",
	"N_B,F_B;inc",
	"N_T,F_T;exc",
	"N_B,F_B;exc",
};

PlotGroup(confs, 45);

//----------------------------------------------------------------------------------------------------
NewRow();

string confs[] = {
	"N_T,F_T;exc",
	"N_B,F_B;exc",
	"N_T,F_T,F_B;exc",
	"N_B,F_T,F_B;exc",
	"N_T,N_B,F_T,F_B;exc",
	"N_T,N_B,F_T;exc",
	"N_T,N_B,F_B;exc",
};

PlotGroup(confs, 25);

GShipout();
