import root;
import pad_layout;

string topDir = "../../";

//string datasets[] = { "DS1", "DS3" };
string datasets[] = { "DS1" };

string arms[] = { "L", "R" };
string units[] = { "N", "F" };
string verts[] = { "T", "B" };

xSizeDef = 8cm;
//ySizeDef = 8cm;


//----------------------------------------------------------------------------------------------------

for (int vi : verts.keys)
{
	NewRow();

	for (int ui : units.keys)
	{
		NewPad("cluster multiplicity (avg.~over planes)", "fraction of events");
		scale(Linear, Log);
		
		int idx = 0;
		for (int dsi : datasets.keys)
		{
			for (int ai : arms.keys)
			{
				string f = topDir+datasets[dsi]+"/showers_combined.root";
				string p = arms[ai]+"/"+units[ui]+"_"+verts[vi]+"/h_avg_mult";
				string l = datasets[dsi] + ", arm " + arms[ai];

				rObject obj = rGetObj(f, p);
				real N = obj.rExec("GetEntries");
				
				draw(shift(0., log10(1./N)), obj, "vl", StdPen(++idx), l);
			}
		}

		limits((-0.5, 0.0001), (30.5, 1), Crop);

		AttachLegend("RP: " + units[ui]+" "+verts[vi]);
	}
}

GShipout();
