import root;
import pad_layout;
include "../run_info.asy";

string topDir = "../../";

string datasets[] = { "DS1" };

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b -- 56t", "45t -- 56b" };
pen dgn_pens[] = { blue, red };

string arms[] = { "_L", "_R" };
string arm_labels[] = { "left arm", "right arm", "double arm" };

xSizeDef = 12cm;
//xTicksDef = LeftTicks(Step=1, step=0.5);


TGraph_errorBar = None;


//----------------------------------------------------------------------------------------------------
NewRow();

for (int ai : arms.keys)
{
	NewRow();

	NewPad("time $\ung{h}$", "mean of $\th_x^*\ung{\mu rad}$");
	//currentpad.yTicks = RightTicks(2., 1);
	DrawRunBands("DS1", -25, +25);
	
	for (int dsi : datasets.keys)
	{
		currentpicture.legend.delete();
		for (int dgi : diagonals.keys)
		{
			string f = topDir+datasets[dsi]+"/distributions_"+diagonals[dgi]+".root";
			pen p = StdPen(dgi+1);
			draw(swToHours*scale(1, 1e6), rGetObj(f, "time dependences/p_th_x"+arms[ai]+"_vs_time"), "eb,d0", p, dgn_labels[dgi]);
		}
	}
	
	limits((2, -25), (3.6, +25), Crop);
	AttachLegend(arm_labels[ai], SE, SE);
	
	for (real y=-25; y <= +25; y += 5)
		xaxis(YEquals(y, false), dotted);
}

GShipout(vSkip=0);
