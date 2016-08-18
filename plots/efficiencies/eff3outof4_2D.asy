import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS1", "DS3" };
//string datasets[] = { "DS1" };

string diagonals[] = { "45b_56t", "45t_56b" };
//string diagonals[] = { "45b_56t" };

string RPs[] = { "L_F", "L_N", "R_N", "R_F" };

string RP_labels[] = { "left far", "left near", "right near", "right far" };


xSizeDef = 6cm;
ySizeDef = 5cm;
//yTicksDef = RightTicks(5., 1.);
//xTicks=LeftTicks(format="$$", 20., 10.)

int gx=0, gy=0;

TH2_palette = Gradient(white, blue, heavygreen, yellow, red);
TH2_z_min = 0;
TH2_z_max = 1.;

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	for (int dgi : diagonals.keys)
	{
		//NewPage();

		string f = topDir + datasets[dsi] + "/eff3outof4_details.root";
		real sgn = (diagonals[dgi] == "45b_56t") ? +1 : -1;
		string opt = "vl,eb";
		
		++gy; gx = 0;
		for (int rpi : RPs.keys)
		{
			++gx;
			NewPad(false, gx, gy);
			label(RP_labels[rpi]);
		}
		
		NewPad(false, -1, gy);
		label(replace("\vbox{\SetFontSizesXX\hbox{dataset: "+datasets[dsi]+"}\hbox{diagonal: "+diagonals[dgi]+"}}", "_", "\_"));

		++gy; gx = 0;
		NewPad(false, -1, gy);	
		label("{\SetFontSizesXII efficiency}");

		for (int rpi : RPs.keys)
		{
			string d = diagonals[dgi] + "/" + RPs[rpi];

			++gx;
			NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$", gx, gy, axesAbove=true);
			draw(scale(1., sgn), rGetObj(f, d+"/track/th_x, th_y : rel"), "def");

			limits((-600, 100), (600, 600), Crop);

			//yaxis(XEquals(-50, false), dashed, above=true);
			//yaxis(XEquals(+80, false), dashed, above=true);
		}
		
		++gy; gx = 0;
		NewPad(false, -1, gy);	
		label("{\SetFontSizesXII more than 1 track}");

		for (int rpi : RPs.keys)
		{
			string d = diagonals[dgi] + "/" + RPs[rpi];

			++gx;
			NewPad("$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$", gx, gy, axesAbove=true);
			draw(scale(1., sgn), rGetObj(f, d+"/pl_suff_no_track/th_x, th_y : rel"), "def");

			limits((-600, 100), (600, 600), Crop);

			//yaxis(XEquals(-50, false), dashed, above=true);
			//yaxis(XEquals(+80, false), dashed, above=true);
		}
	}
}

GShipout(vSkip=0pt);
