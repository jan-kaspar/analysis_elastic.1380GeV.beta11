import root;
import pad_layout;

string topDir = "../../";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

string file_45b = topDir + "DS1/distributions_45b_56t.root";
string file_45t = topDir + "DS1/distributions_45t_56b.root";

string rps[] = { "L_F", "L_N", "R_N", "R_F" };
string rp_labels[] = { "left far", "left near", "right near", "right far" };

real sh_top[] = { +3.60, +3.20, +2.70, +2.90 };
real sh_bot[] = { -2.95, -2.75, -2.90, -3.10 };

//----------------------------------------------------------------------------------------------------

real edge = 3.6101;
real cutEdge = 2.22721 / sqrt(2);
int strips = 11;
real margin_v_e = 0.2;
real margin_v_b = 0.4;
real margin_u = 0.1;

path det_shape = (cutEdge, 0)--(edge, 0)--(edge, edge)--(0, edge)--(0, cutEdge)--cycle;
det_shape = scale(10) * rotate(45) * shift(-cutEdge/2, -cutEdge/2) * det_shape;

path hor_det_shape = shift(0, -cutEdge/sqrt(2)*10) * det_shape;

//----------------------------------------------------------------------------------------------------

xTicksDef = LeftTicks(Step=10, step=5);
yTicksDef = RightTicks(Step=10, step=5);
//yTicksDef = RightTicks(Step=0.5, step=0.1);

NewPad(false);
label("{\SetFontSizesXX before selection}");

for (int ri : rps.keys)
{
	NewPad("$x\ung{mm}$", "$y\ung{mm}$", 6cm, 10cm);
	scale(Linear, Linear, Log);

	TH2_z_max = log10(2e3);
	
	draw(rGetObj(file_45b, "hit distributions/vertical, aligned, before selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_nosel"), "p,bar");
	draw(rGetObj(file_45t, "hit distributions/vertical, aligned, before selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_nosel"), "p");
	
	draw(shift(0, sh_top[ri])*det_shape);
	draw(shift(0, sh_bot[ri])*scale(1, -1)*det_shape);

	limits((-30, -50), (+30, +50), Crop);
	//limits((-30, -5), (+30, +5), Crop);

	for (real x = -30; x <= +30; x += 5)
		yaxis(XEquals(x, false), dotted+gray);

	for (real y = -50; y <= +50; y += 5)
		xaxis(YEquals(y, false), dotted+gray);
	//for (real y = -5; y <= +5; y += 0.2)
	//	xaxis(YEquals(y, false), dotted+gray);

	AttachLegend(replace(rp_labels[ri], "_", "\_"));
}

//----------------------------------------------------------------------------------------------------

NewRow();

NewPad(false);
label("{\SetFontSizesXX after selection}");

for (int ri : rps.keys) {
	NewPad("$x\ung{mm}$", "$y\ung{mm}$", 6cm, 10cm);
	scale(Linear, Linear, Log);

	TH2_z_max = log10(6e2);
	
	draw(rGetObj(file_45b, "hit distributions/vertical, aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_sel"), "p,bar");
	draw(rGetObj(file_45t, "hit distributions/vertical, aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_sel"), "p");
	
	draw(shift(0, sh_top[ri])*det_shape);
	draw(shift(0, sh_bot[ri])*scale(1, -1)*det_shape);

	limits((-30, -50), (+30, +50), Crop);
	//limits((-30, -5), (+30, +5), Crop);

	for (real x = -30; x <= +30; x += 5)
		yaxis(XEquals(x, false), dotted+gray);

	for (real y = -50; y <= +50; y += 5)
		xaxis(YEquals(y, false), dotted+gray);

	AttachLegend(replace(rp_labels[ri], "_", "\_"));
}
