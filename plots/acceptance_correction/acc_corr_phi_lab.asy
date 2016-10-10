import root;
import pad_layout;

string topDir = "../../";

string dataSets[] = {
	"DS1",
//	"DS3"
};

real th_y_low_45b[] = { 210, 500 };
real th_y_low_45t[] = { -210, -500 };

real radii[] = { 200, 300, 400, 500 };

TH2_palette = Gradient(blue, heavygreen, yellow, red);

for (int dsi : dataSets.keys)
{
	NewPad("$\th_x^{*}\ung{\mu rad}$", "$\th_y^{*}\ung{\mu rad}$");
	scale(Linear, Linear, Log);
	
	label("$\th^* (\mu rad)$:", (radii[0], 0), W);
	//label(rotate(-90)*Label("\SmallerFonts $\rm\mu rad$"), (165, 0), 0.5E);
	
	for (int ri : radii.keys)
	{
		draw(scale(radii[ri])*unitcircle, dashed);
		label(rotate(-90)*Label("\SmallerFonts $"+format("%.0f", radii[ri])+"$"), (radii[ri], 0), 0.5E);
	}

	
	// 45 bottom - 56 top
	TH2_z_min = 10^5;
	TH2_z_max = 10^7;
	draw(scale(1e6, 1e6), RootGetObject(topDir + dataSets[dsi]+"/distributions_45b_56t.root", "normalization/h_th_y_vs_th_x_normalized"), "def");
	
	real th_y_low = th_y_low_45b[dsi], th_y_high = 600;
	
	draw((-600, th_y_low)--(600, th_y_low), magenta+1pt);
	//draw((-600, th_y_high)--(600, th_y_high), magenta+1pt);
	
	for (int ri : radii.keys)
	{
		if (radii[ri] > th_y_low)
			draw(arc((0, 0), radii[ri], aSin(th_y_low/radii[ri]), 180 - aSin(th_y_low/radii[ri])), black+1pt);
	}
	
	// 45 top - 56 bottom
	//TH2_z_min = 7;
	//TH2_z_max = 9.6;
	draw(scale(1e6, 1e6), RootGetObject(topDir + dataSets[dsi]+"/distributions_45t_56b.root", "normalization/h_th_y_vs_th_x_normalized"), "p");
	
	real th_y_low = th_y_low_45t[dsi], th_y_high = -600;
	
	draw((-600, th_y_low)--(600, th_y_low), magenta+1pt);
	//draw((-600, th_y_high)--(600, th_y_high), magenta+1pt);
	
	for (int ri : radii.keys)
	{
		if (radii[ri] > abs(th_y_low))
			draw(arc((0, 0), radii[ri], 180 - aSin(th_y_low/radii[ri]), 360 + aSin(th_y_low/radii[ri])), black+1pt);
	}
	
	limits((-600, -600), (600, 600), Crop);
	AttachLegend(dataSets[dsi]);
}
