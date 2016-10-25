import root;
import pad_layout;

string rps[] = {
	"L_F",
	"L_N",
	"R_N",
	"R_F",
};

string f = "test1.root";

TGraph_errorBar = None;

xSizeDef = 8cm;

xTicksDef = LeftTicks(200., 100.);

//----------------------------------------------------------------------------------------------------

frame f_leg;

for (int rpi : rps.keys)
{
	NewPad("$\th^*_{x,y}\ung{\mu rad}$");
	scale(Linear, Log);

	real th_x_sh = 1.;

	draw(scale(1e6, 1)*shift(0, th_x_sh), RootGetObject(f, "th_x/45b_56t"), "p", mCi+1pt+black, "$\th_x^*$, histogram");
	draw(scale(1e6, 1)*shift(0, th_x_sh), RootGetObject(f, "th_x/45b_56t|ff"), red+1pt, "$\th_x^*$, fit");

	draw(scale(1e6, 1), RootGetObject(f, "th_y/"+rps[rpi]), "p", heavygreen, mCi+1pt+heavygreen, "$\th_y^*$, histogram");
	draw(scale(1e6, 1), RootGetObject(f, "th_y/"+rps[rpi]+"|ff"), blue+1pt, "$\th_y^*$, fit");

	limits((-600., 10.^-1), (+600, 10.^4), Crop);

	f_leg = BuildLegend();

	currentpicture.legend.delete();

	AttachLegend(replace(rps[rpi], "_", "\_"));
}

NewPad(false);
add(f_leg);
