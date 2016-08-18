import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS1", "DS3" };
//string datasets[] = { "DS1" };

string diagonals[] = { "45b_56t", "45t_56b" };
//string diagonals[] = { "45b_56t" };

string RPs[] = { "L_F", "L_N", "R_N", "R_F" };

string RP_labels[] = { "left far", "left near", "right near", "right far" };


xSizeDef = 8cm;
ySizeDef = 8cm;
yTicksDef = RightTicks(10., 2.);

//legendLabelPen = fontcommand("\SetFontSizesXX");

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	NewPage();
	int gx=0, gy=0;

	for (int dgi : diagonals.keys)
	{
		string f = topDir + datasets[dsi] + "/eff3outof4_details.root";
		real sgn = (diagonals[dgi] == "45b_56t") ? +1 : -1;

		string opt = "d0,eb";
		
		++gy; gx = 0;
		for (int rpi : RPs.keys)
		{
			++gx;
			NewPad(false, gx, gy);
			label("{\SetFontSizesXX "+RP_labels[rpi]+"}");
		}
		
		++gy;
		NewPad(false, -1, gy);
		label(replace("\vbox{\SetFontSizesXX\hbox{dataset: "+datasets[dsi]+"}\hbox{diagonal: "+diagonals[dgi]+"}}", "_", "\_"));

		gx = 0;
		for (int rpi : RPs.keys)
		{
			string d = diagonals[dgi] + "/" + RPs[rpi];

			real w = 1.2pt;

			++gx;
			NewPad("$\th_y^*\ung{\mu rad}$", "\ung{\%}", gx, gy);
			//currentpad.xTicks=LeftTicks(format="$$", 20., 10.);
			draw(scale(sgn, 100), rGetObj(f, d+"/anything/th_y : rel"), opt, black+w);
			draw(scale(sgn, 100), rGetObj(f, d+"/track/th_y : rel"), opt, cyan+w);
			draw(scale(sgn, 100), rGetObj(f, d+"/track_compatible/th_y : rel"), opt, magenta+w);
			
			draw(scale(sgn, 100), rGetObj(f, d+"/pl_insuff/th_y : rel", error=false), opt, red+w);
			draw(scale(sgn, 100), rGetObj(f, d+"/pl_suff_no_track/th_y : rel"), opt, blue+w);
			draw(scale(sgn, 100), rGetObj(f, d+"/pat_more/th_y : rel"), opt, heavygreen+w);

			limits((100, 0), (600, 100), Crop);

			for (real y = 0; y <= 100; y += 10)
				xaxis(YEquals(y, false), dotted);
		}
		
		/*
		++gy; gx = 0;
		for (int rpi : RPs.keys) {
			string d = diagonals[dgi] + "/" + RPs[rpi];
			
			++gx;
			NewPad("$\th_y^*\ung{\mu rad}$", "\ung{\%}", gx, gy);
			draw(scale(sgn, 100), rGetObj(f, d+"/pl_insuff/th_y : rel", error=false), opt, red);
			draw(scale(sgn, 100), rGetObj(f, d+"/pl_suff_no_track/th_y : rel"), opt, blue);
			draw(scale(sgn, 100), rGetObj(f, d+"/pat_more/th_y : rel"), opt, heavygreen);

			limits((0, 0), (600, 40), Crop);
		}
		*/

		NewPad(false, gx+1, gy);
		AddToLegend("anything", black);
		AddToLegend("track", cyan);
		AddToLegend("compatible track", magenta);
		AddToLegend("pl\_insuff", red);
		AddToLegend("pl\_suff\_no\_track", blue);
		AddToLegend("pat\_more", heavygreen);
		AttachLegend();
	}
}

GShipout(vSkip=0pt, hSkip=1mm);
