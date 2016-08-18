import root;
import pad_layout;

string topDir = "../../";

string f = topDir + "DS1/alignment.root";
string dir = "period 1/unit L_F";

transform xyswitch = (0, 0, 0, 1, 1, 0);

NewPad("$y\ung{mm}$", "$x\ung{mm}$", 12cm, 6cm);
draw(rGetObj(f, dir+"/horizontal/horizontal graph fit/horizontal fit|merged"), "p", heavygreen);
draw(rGetObj(f, dir+"/horizontal/horizontal profile/p"), "eb,d0", red+1pt);
draw(rGetObj(f, dir+"/horizontal/horizontal profile/p|ff"), blue+1pt);

label("bottom RP", (-6, -1.5));
label("top RP", (+6, -1.5));

limits((-12, -2), (+12, +2), Crop);
xaxis(YEquals(0, false), dashed);

GShipout("alignment_method_x");

//----------------------------------------------------------------------------------------------------

string f = topDir + "DS1/alignment.root";
string dir = "period 2/unit L_F";


NewPad("$y\ung{mm}$", "", 12cm, 6cm);
draw(rGetObj(f, dir+"/vertical/y_hist"), "vl", red);

label("bottom RP", (-8, 500));
label("top RP", (+8, 500));

limits((-12, 0), (+12, 900), Crop);

GShipout("alignment_method_y");
