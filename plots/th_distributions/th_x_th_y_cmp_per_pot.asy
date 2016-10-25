import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	//"DS1_no_add_alignment",
	"DS1",
	//"DS1_glob_al_45t_56b",
};

string dgns[] = {
	"45b_56t",
	"45t_56b"
};

string rps[] = {
	"L_F",
	"L_N",
	"R_N",
	"R_F",
};

xSizeDef = 10cm;

xTicksDef = LeftTicks(100., 50.);

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	NewRow();

	frame f_leg;

	for (int rpi : rps.keys)
	{
		if (rpi == 2)
			NewRow();

		NewPad("$\th_{x,y}^*\ung{\mu rad}$");
		scale(Linear, Log);
	
		int ci = -1;
	
		for (int dgni : dgns.keys)
		{
			string dgn_label = replace(dgns[dgni], "_", " -- ");
	
			transform tr = scale(1e6, 1);
	
			string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgni] + ".root";
	
			string opt = "vl";
			//string opt = "eb";
	
			TH1_x_min = -inf; TH1_x_max = +inf;
			draw(shift(0, 1.15)*tr, RootGetObject(f, "selected - angles/h_th_x"), opt, StdPen(++ci), dgn_label + ", $\th_x^*$ (scaled)");

			if (dgns[dgni] == "45b_56t")
				TH1_x_min = +230e-6; 
			else
				TH1_x_max = -230e-6;
			draw(tr, RootGetObject(f, "selected - angles/h_th_y_" + rps[rpi]), opt, StdPen(++ci)+1.5pt, dgn_label + ": $\th_y^*$");
		}
	
		f_leg = BuildLegend(replace(datasets[dsi], "_", "\_"));

		currentpicture.legend.delete();

		limits((-500, 0.1), (+500, 1e4), Crop);
		AttachLegend(replace(rps[rpi], "_", "\_"));
	}

	NewPad(false);
	add(f_leg);
}
