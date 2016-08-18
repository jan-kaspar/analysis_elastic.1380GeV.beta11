import pad_layout;

real si_vtx_y = 70e-6;

void DrawCut(bool leftArm, real v_y, real L_y, real y, pen p=black, string label)
{
	real vtx_y_min = -4 * si_vtx_y;
	real vtx_y_max = +4 * si_vtx_y;
	
	real th_y_min, th_y_max;
	if (leftArm)
	{
		th_y_min = (-y + v_y * vtx_y_min) / L_y;
		th_y_max = (-y + v_y * vtx_y_max) / L_y;
	} else {
		th_y_min = (y - v_y * vtx_y_min) / L_y;
		th_y_max = (y - v_y * vtx_y_max) / L_y;
	}

	draw(scale(1e6, 1e6) * ((vtx_y_min, th_y_min)--(vtx_y_max, th_y_max)), p, label);
}

//----------------------------------------------------------------------------------------------------

// 45 bot -- 56 top
NewPad("$y^*\ung{\mu m}$", "$\th_y^{*L,R}\ung{\mu rad}$");
currentpad.yTicks = RightTicks(50., 10.);
DrawCut(true, -3.3396, 21.02417, -2.95e-3, black, "L, F");
DrawCut(true, -3.0864, 21.03405, -2.75e-3, red, "L, N");
DrawCut(false, -3.0412, 18.46229, +2.70e-3, blue, "R, N");
DrawCut(false, -3.2815, 18.15756, +2.90e-3, heavygreen, "R, F");
AttachLegend("45 bot -- 56 top", NW, NE);

GShipout(margin=1mm);
