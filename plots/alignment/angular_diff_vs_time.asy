import root;
import pad_layout;

string topDir = "../../";

//string datasets[] = { "DS1", "DS3" };
string datasets[] = { "DS1" };

string dgns[] = { "45b_56t", "45t_56b" };
string dgn_labs[] = { "45 bot -- 56 top", "45 top -- 56 bot" };


xSizeDef = 18cm;
xTicksDef = LeftTicks(Step = 0.5, step = 0.1);

real time_min = 2, time_max = 3.6;

//----------------------------------------------------------------------------------------------------

NewPad(false);
label("{\SetFontSizesXX F-N difference, left arm}");

NewPad("time $\ung{h}$", "$\De^{F-N} \th_y^{*L}\ung{\mu rad}$");
for (int dsi : datasets.keys)
{
	for (int dgi : dgns.keys)
	{
		string f = topDir+datasets[dsi]+"/distributions_"+dgns[dgi]+".root";
		draw(scale(1/3600, 1e6), RootGetObject(f, "time dependences/p_diffNF_th_y_L_vs_time"), "eb,d0", StdPen(dgi+1), dgn_labs[dgi]);
	}
}
limits((time_min, -10), (time_max, +2), Crop);
xaxis(YEquals(0, false), dotted);
AttachLegend(NW, NE);

//----------------------------------------------------------------------------------------------------
NewRow();

NewPad(false);
label("{\SetFontSizesXX F-N difference, right arm}");

NewPad("time $\ung{h}$", "$\De^{F-N} \th_y^{*R}\ung{\mu rad}$");
for (int dsi : datasets.keys)
{
	for (int dgi : dgns.keys)
	{
		string f = topDir+datasets[dsi]+"/distributions_"+dgns[dgi]+".root";
		draw(scale(1/3600, 1e6), RootGetObject(f, "time dependences/p_diffNF_th_y_R_vs_time"), "eb,d0", StdPen(dgi+1));
	}
}

limits((time_min, -0.2), (time_max, +1.1), Crop);
xaxis(YEquals(0, false), dotted);

//----------------------------------------------------------------------------------------------------
NewRow();

NewPad(false);
label("{\SetFontSizesXX R-L difference, extraploated}");

TGraph_errorBar = None;
NewPad("time $\ung{h}$", "$\De^{R-L} \th_y^{*}\ung{\mu rad}$");
for (int dsi : datasets.keys)
{
	for (int dgi : dgns.keys)
	{
		string f = topDir+datasets[dsi]+"/distributions_"+dgns[dgi]+".root";
		draw(scale(1/3600, 1e6), RootGetObject(f, "time dependences/g_ext_diffLR_th_y_vs_time"), "p", StdPen(dgi+1), mCi+2pt+StdPen(dgi+1));
	}
}

limits((time_min, -20), (time_max, +30), Crop);
xaxis(YEquals(0, false), dotted);


//----------------------------------------------------------------------------------------------------
NewRow();

NewPad(false);
label("{\SetFontSizesXX R-L difference, extraploated}");

TGraph_errorBar = None;
NewPad("time $\ung{h}$", "$\De^{R-L} \th_x^{*}\ung{\mu rad}$");
for (int dsi : datasets.keys)
{
	for (int dgi : dgns.keys)
	{
		string f = topDir+datasets[dsi]+"/distributions_"+dgns[dgi]+".root";
		draw(scale(1/3600, 1e6), RootGetObject(f, "time dependences/g_ext_diffLR_th_x_vs_time"), "p", StdPen(dgi+1), mCi+2pt+StdPen(dgi+1));
	}
}

limits((time_min, -2), (time_max, +3), Crop);
xaxis(YEquals(0, false), dotted);
