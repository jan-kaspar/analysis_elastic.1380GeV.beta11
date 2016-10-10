import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	"DS1_no_add_alignment",
	"DS1",
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
		NewPad("$\th_{x,y}^*\ung{\mu rad}$");
		scale(Linear, Log);
	
		int ci = -1;
	
		for (int dgni : dgns.keys)
		{
			string dgn_label = replace(dgns[dgni], "_", " -- ");
	
			transform tr = scale(1e6, 1);
	
			string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgni] + ".root";
	
			//string opt = "vl";
			string opt = "eb";
	
			draw(shift(0, 1.3)*tr, RootGetObject(f, "selected - angles/h_th_x"), opt, StdPen(++ci), dgn_label + ", $\th_x^*$");

			draw(tr, RootGetObject(f, "selected - angles/h_th_y_" + rps[rpi]), opt, StdPen(++ci), dgn_label + ": $\th_y^*$");
		}
	
		f_leg = BuildLegend(replace(datasets[dsi], "_", "\_"));

		currentpicture.legend.delete();

		limits((-500, 0.1), (+500, 1e4), Crop);
		AttachLegend(replace(rps[rpi], "_", "\_"));
	}

	NewPad(false);
	add(f_leg);
}
