import root;
import pad_layout;
include "../run_info.asy";

string topDir = "../../";

string dataset = "DS1";

string units[] = { "L_F", "L_N", "R_N", "R_F" };
string unit_labels[] = { "left, far", "left, near", "right, near", "right, far" };

/*
string units[] = { "L_1_F", "L_1_N", "R_1_F" };
string unit_labels[] = { "left, 210, far", "left, 210, near", "right, 210, far" };
*/

xSizeDef = 10cm;
drawGridDef = true;

TGraph_errorBar = None;

string period = "period 3";

//----------------------------------------------------------------------------------------------------
NewRow();

for (int ui : units.keys)
{
	NewPad("$y\ung{mm}$", "$\hbox{mean } x\ung{mm}$");
	currentpad.xTicks = LeftTicks(2., 1.);

	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/horizontal/horizontal profile/p"), "d0,eb", blue);
	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/horizontal/horizontal profile/p|ff"), "l", red+1pt);
	
	limits((-10, -1), (+10, +1), Crop);
	AttachLegend(unit_labels[ui], NE, NE);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int ui : units.keys)
{
	NewPad("$y\ung{mm}$", "");
	currentpad.xTicks = LeftTicks(2., 1.);
	scale(Linear, Log);

	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/vertical/y_hist|y_hist"), "d0,vl", blue);
	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/vertical/y_hist|y_hist_range"), "d0,vl", red);
	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/vertical/y_hist|y_hist_range|gaus"), "d0,vl", heavygreen+dashed+1pt);

	limits((-10, 1), (+10, 1e3), Crop);
	AttachLegend(unit_labels[ui], NE, NE);
	xaxis(YEquals(3, false), black+1pt);
}


//----------------------------------------------------------------------------------------------------
NewRow();

for (int ui : units.keys)
{
	NewPad("$\de y\ung{mm}$", "");

	if (units[ui] == "R_1_N")
		continue;

	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/vertical/g_max_diff"), "l,p", heavygreen, mCi+1pt+heavygreen);

	limits((-2.5, 0), (+2.5, 1), Crop);
	AttachLegend(unit_labels[ui], NE, NE);
}
