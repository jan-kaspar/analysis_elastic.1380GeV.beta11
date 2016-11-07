import root;
import pad_layout;

string topDir = "./";

string dir_45b = "45b_56t/L_F/sh = -0.00";
string dir_45t = "45t_56b/L_F/sh = -20.00";


xSizeDef = 8cm;

xTicksDef = LeftTicks(100., 50.);

//----------------------------------------------------------------------------------------------------

string f = topDir + "match_th_y_scale_fcn.root";

NewPad("$\th^*_{x,y}\ung{\mu rad}$", "");
scale(Linear, Log);

real sc_45b = 1.;
//real sc_45t = 1.025;
real sc_45t = 1.;

draw(scale(1e6, sc_45b), RootGetObject(f, dir_45b+"/h_cmp|g_h_th_x"), "l", black);
draw(scale(1e6, sc_45b), RootGetObject(f, dir_45b+"/h_cmp|h_th_y"), "vl", red);

draw(scale(-1e6, sc_45t), RootGetObject(f, dir_45t+"/h_cmp|g_h_th_x"), "l", blue);
draw(scale(-1e6, sc_45t), RootGetObject(f, dir_45t+"/h_cmp|h_th_y"), "vl", magenta);
