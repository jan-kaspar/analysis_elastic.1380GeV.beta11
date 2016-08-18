import root;
import pad_layout;
include "../run_info.asy";

string topDir = "../../";

//string datasets[] = { "DS1", "DS3" };
string datasets[] = { "DS1" };
real timestamp0[] = { 1360537200, 1360796400 };

string units[] = { "L_F", "L_N", "R_N", "R_F" };
string bpms[] = { "BPMWT.B6L5.B2", "BPMWT.A6L5.B2", "BPMWT.A6R5.B1", "BPMWT.B6R5.B1" };
real sh_x[] = { 250, -250, 870, 400 };
real sh_y[] = { +200, +660, -300, 0 };

/*
string units[] = { "L_N", "R_N" };
string bpms[] = { "BPMWT.A6L5.B2", "BPMWT.A6R5.B1" };
real sh_x[] = { -350, +770 };
real sh_y[] = { +750, -200 };
*/

xSizeDef = 8cm;
ySizeDef = 6cm;

drawGridDef = true;

TGraph_errorBar = None;

real t_min = 1.5, t_max = 4.0;


NewPad(false);
label("\SetFontSizesXX\vbox{\hbox{horizontal}\hbox{shift}}");

for (int ui : units.keys)
{
	NewPad("time $\ung{h}$", "horizontal position $\ung{\mu m}$", yTicks=RightTicks(Step=50, step=10));
	DrawRunBands(-100, +100);

	/*
	TGraph_reducePoints = 30;
	draw(shift(-timestamp0, sh_x[ui]), RootGetObject("bpm.root", "LHC.BOFSU:POSITIONS_H::"+bpms[ui]), black);
	TGraph_reducePoints = 1; 
	*/

	for (int di : datasets.keys)
	{
		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/b_p"), "p,l,eb", cyan, mCi+1pt+cyan);
		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/b_g"), "p,l,eb", green, mCi+1pt+green);

		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/b"), "p,l,eb", blue+1pt, mCi+1pt+blue);

		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment_fit.root", ""+units[ui]+"/b_fit"), "l", red+1.5pt);
	}

	limits((t_min, -100), (t_max, +100), Crop);
	AttachLegend(replace(units[ui], "_", "\_"), SE, SE);
}

//----------------------------------------------------------------------------------------------------
//NewPage();
NewRow();

NewPad(false);
label("\SetFontSizesXX\vbox{\hbox{vertical}\hbox{shift}}");

for (int ui : units.keys)
{
	NewPad("time $\ung{h}$", "vertical position $\ung{\mu m}$");
	DrawRunBands(-200, +500);

	//TGraph_reducePoints = 30;
	//draw(shift(-timestamp0, sh_y[ui]), RootGetObject("bpm.root", "LHC.BOFSU:POSITIONS_V::"+bpms[ui]), black);
	//TGraph_reducePoints = 1; 

	for (int di : datasets.keys)
	{
		pen p = StdPen(di+1);
		//draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/c_min_diff"), "p,l,eb", cyan, mCi+1pt+cyan);
		//draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/c_prob"), "p,l,eb", green, mCi+1pt+green);
		//draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/c_mean_diff_sq"), "p,l,eb", magenta, mCi+1pt+magenta);
		//draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/c_hist_chi_sq"), "p,l,eb", green, mCi+1pt+green);
		
		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/c"), "p,l,eb", blue+1pt, mCi+1pt+blue);

		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment_fit.root", ""+units[ui]+"/c_fit"), "l", red+1.5pt);
	}

	limits((t_min, -200), (t_max, +500), Crop);
	AttachLegend(replace(units[ui], "_", "\_"), SE, SE);
}

//----------------------------------------------------------------------------------------------------
NewRow();

NewPad(false);
label("\SetFontSizesXX\vbox{\hbox{tilt in}\hbox{$xy$ plane}}");

for (int ui : units.keys)
{
	NewPad("time $\ung{h}$", "tilt $\ung{mrad}$", yTicks=RightTicks(Step=5, step=1));
	DrawRunBands(-5, +5);

	for (int di : datasets.keys)
	{
		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/a_p"), "p,l,eb", cyan, mCi+1pt+cyan);
		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/a_g"), "p,l,eb", green, mCi+1pt+green);
		
		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment.root", "global/"+units[ui]+"/a"), "p,l,eb", blue, mCi+1pt+blue);

		draw(swToHours, RootGetObject(topDir + datasets[di]+"/alignment_fit.root", ""+units[ui]+"/a_fit"), "l", red+1.5pt);
	}

	limits((t_min, -5), (t_max, +5), Crop);
	AttachLegend(replace(units[ui], "_", "\_"), SE, SE);
}

GShipout();
