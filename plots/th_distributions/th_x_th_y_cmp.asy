import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS1" };

string dgns[] = {
	"45b_56t",
	"45t_56b"
};

xSizeDef = 10cm;

xTicksDef = LeftTicks(100., 50.);

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	NewPad("$\th_{x,y}^*\ung{\mu rad}$");
	scale(Linear, Log);

	int ci = -1;

	for (int dgni : dgns.keys)
	{
		string dgn_label = replace(dgns[dgni], "_", " -- ");

		transform tr = scale(1e6, 1);

		string f = topDir + datasets[dsi] + "/distributions_" + dgns[dgni] + ".root";

		string opt = "vl";

		draw(shift(0, 1.3)*tr, RootGetObject(f, "selected - angles/h_th_x"), opt, StdPen(++ci), dgn_label + ", $\th_x^*$");
		draw(tr, RootGetObject(f, "selected - angles/h_th_y"), opt, StdPen(++ci), dgn_label + ": $\th_y^*$");
	}

	AttachLegend(NW, NE);
}
