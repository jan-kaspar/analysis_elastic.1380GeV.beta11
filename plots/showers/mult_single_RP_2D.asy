import root;
import pad_layout;

string topDir = "../../";

//string datasets[] = { "DS1", "DS3" };
string datasets[] = { "DS1" };

string arms[] = { "L", "R" };
string units[] = { "N", "F" };
string verts[] = { "T", "B" };

string N_pl = "3";

xSizeDef = 8cm;
ySizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

NewPad(false);

for (int dsi : datasets.keys)
{
	for (int ai : arms.keys)
	{
		NewPad(false);
		label("{\SetFontSizesXX " + datasets[dsi] + ", arm " + arms[ai] + "}");
	}
}

for (int ui : units.keys)
{
	for (int vi : verts.keys)
	{
		NewRow();

		NewPad(false);
		label("{\SetFontSizesXX " + units[ui]+" "+verts[vi] + "}");
		
		for (int dsi : datasets.keys)
		{
			for (int ai : arms.keys)
			{
				NewPad("avg.~multiplicity over first " + N_pl + " planes", "avg.~multiplicity over last " + N_pl + " planes");
				scale(Linear, Linear, Log);

				string f = topDir+datasets[dsi]+"/showers_combined.root";
				string p = arms[ai]+"/"+units[ui]+"_"+verts[vi]+"/h2_avg_mult_fl_pl_" + N_pl;
				draw(rGetObj(f, p));

				draw((0, 0)--(30, 30), magenta+2pt);
			}
		}
	}
}

GShipout();
