import root;
import pad_layout;

string topDir = "../../";

string f = topDir + "reconstruction_formulae/test_formulae_graph.root";

void DrawSet(string caption, string iq, real scale_x, string unit_x, string unit_y, string formulae[])
{
	NewPad(false, 0, -1);
	label("{\SetFontSizesXX " + caption + "}");

	string complementary = (iq == "x" || iq == "y") ? "angle" : "vertex";

	for (int fi : formulae.keys)
	{
		string formula = formulae[fi];
	
		NewRow();
	
		NewPad(false);
		label("{\SetFontSizesXX " + replace(formula, "_", "\_") + "}");
		
		NewPad("$"+iq+"^{*, \rm sim}\ung{"+unit_x+"}$", "std.~dev.~of $"+iq+"^{*, \rm reco} - "+iq+"^{*, \rm sim}\ung{"+unit_y+"}$");
		scale(Linear, Linear(true));
		draw(scale(scale_x, 1), rGetObj(f, formula+"/pitch/g_stddev"), "l,p", black, mCi+1pt+black, "pitch");
		draw(scale(scale_x, 1), rGetObj(f, formula+"/beamDiv/g_stddev"), "l,p", red, mCi+1pt+red, "beam div.");
		draw(scale(scale_x, 1), rGetObj(f, formula+"/"+complementary+"/g_stddev"), "l,p", blue, mCi+1pt+blue, ""+complementary+"");
		draw(scale(scale_x, 1), rGetObj(f, formula+"/pitch,beamDiv,"+complementary+"/g_stddev"), "l,p", magenta, mCi+1pt+magenta, "pitch, beam div., "+complementary+"");

		if (fi == 0)
		{
			frame f_leg = BuildLegend();
			NewPad(false, 1, -1);
			attach(f_leg);	
		}
		
		NewPad("$"+iq+"^{*, \rm sim}\ung{"+unit_x+"}$", "std.~dev.~of $"+iq+"^{*, \rm reco} - "+iq+"^{*, \rm sim}\ung{"+unit_y+"}$");
		scale(Linear, Linear(true));
		draw(scale(scale_x, 1), rGetObj(f, formula+"/optics/g_stddev"), "l,p", heavygreen, mCi+1pt+heavygreen, "optics");
		draw(scale(scale_x, 1), rGetObj(f, formula+"/"+complementary+",optics/g_stddev"), "l,p", brown, mCi+1pt+brown, ""+complementary+", optics");

		if (fi == 0)
		{
			frame f_leg = BuildLegend();
			NewPad(false, 2, -1);
			attach(f_leg);	
		}
		
		NewPad("$"+iq+"^{*, \rm sim}\ung{"+unit_x+"}$", "std.~dev.~of $"+iq+"^{*, \rm reco} - "+iq+"^{*, \rm sim}\ung{"+unit_y+"}$");
		scale(Linear, Linear(true));
		draw(scale(scale_x, 1), rGetObj(f, formula+"/pitch,beamDiv,"+complementary+",optics/g_stddev"), "l,p", orange+1pt, mCi+1.5pt+orange, "pitch, beam div., "+complementary+", optics");

		if (fi == 0)
		{
			frame f_leg = BuildLegend();
			NewPad(false, 3, -1);
			attach(f_leg);	
		}
	}
}

//----------------------------------------------------------------------------------------------------

string formulae[] = {
	"theta_x,1_arm,hit:th_x_L",
	"theta_x,1_arm,angle:th_x_L",
	"theta_x,1_arm,regr:th_x_L",
};

DrawSet("$\th_x^*$, single arm", "\th_x", 1e6, "\mu rad", "\mu rad", formulae);

//----------------------------------------------------------------------------------------------------
NewPage();

string formulae[] = {
	"theta_x,2_arm,LRavg_hit:th_x",
	"theta_x,2_arm,LRavg_angle:th_x",
	"theta_x,2_arm,LRavg_regr:th_x",
	"theta_x,2_arm,regr:th_x",
};

DrawSet("$\th_x^*$, double arm", "\th_x", 1e6, "\mu rad", "\mu rad", formulae);

//----------------------------------------------------------------------------------------------------
NewPage();

string formulae[] = {
	"theta_y,1_arm,one_pot_hit_LF:th_y_L",
	"theta_y,1_arm,one_pot_hit_LN:th_y_L",
};

DrawSet("$\th_y^*$, single arm, single RP", "\th_y", 1e6, "\mu rad", "\mu rad", formulae);

//----------------------------------------------------------------------------------------------------
NewPage();

string formulae[] = {
	"theta_y,1_arm,hit:th_y_L",
	"theta_y,1_arm,angle:th_y_L",
	"theta_y,1_arm,regr:th_y_L",
};

DrawSet("$\th_y^*$, single arm", "\th_y", 1e6, "\mu rad", "\mu rad", formulae);

//----------------------------------------------------------------------------------------------------
NewPage();
	
string formulae[] = {
	"theta_y,2_arm,LRavg_hit:th_y",
	"theta_y,2_arm,LRavg_angle:th_y",
	"theta_y,2_arm,LRavg_regr:th_y",
	"theta_y,2_arm,regr:th_y",
};

DrawSet("$\th_y^*$, double arm", "\th_y", 1e6, "\mu rad", "\mu rad", formulae);

//----------------------------------------------------------------------------------------------------
NewPage();

string formulae[] = {
	"vtx_x,1_arm,regr:vtx_x_L",
};

DrawSet("$x^*$, single arm", "x", 1e3, "\mu m", "\mu m", formulae);

//----------------------------------------------------------------------------------------------------
NewPage();

string formulae[] = {
	"vtx_x,2_arm,LRavg_regr:vtx_x",
	"vtx_x,2_arm,regr:vtx_x",
};

DrawSet("$x^*$, double arm", "x", 1e3, "\mu m", "\mu m", formulae);

//----------------------------------------------------------------------------------------------------
NewPage();

string formulae[] = {
	"vtx_y,1_arm,regr:vtx_y_L",
};

DrawSet("$y^*$, single arm", "y", 1e3, "\mu m", "\mu m", formulae);

//----------------------------------------------------------------------------------------------------
NewPage();

string formulae[] = {
	"vtx_y,2_arm,LRavg_regr:vtx_y",
	"vtx_y,2_arm,regr:vtx_y",
};

DrawSet("$y^*$, double arm", "y", 1e3, "\mu m", "\mu m", formulae);
