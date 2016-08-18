import root;
import pad_layout;

string topDir = "../../";

//string datasets[] = { "DS1", "DS3" };
string datasets[] = { "DS1" };

//string diagonals[] = { "45b_56t", "45t_56b" };
string diagonals[] = { "45b_56t" };

string elements[][] = {
	{ "L_N/#", "L_F/#", "L_N, L_F/# && #" },
	{ "R_N/#", "R_F/#", "R_N, R_F/# && #" },
	{ "", "", "dgn/# && #, L || R" }
};


//----------------------------------------------------------------------------------------------------

TGraph_errorBar = None;


for (int dsi : datasets.keys)
{

	for (int di : diagonals.keys)
	{
		string dgn = diagonals[di];

		string f = topDir + datasets[dsi]+"/distributions_"+dgn+".root";
		
		NewPage();

		NewPad("$\th_y^{*,LN} \ung{\mu rad}$", "$\th_y^{*,LF} \ung{\mu rad}$");
		draw(scale(1e6, 1e6), RootGetObject(f, "optics/p_th_y_LF_vs_th_y_LN"), "d0,eb", red);
		draw(scale(1e6, 1e6), RootGetObject(f, "optics/p_th_y_LF_vs_th_y_LN|pol1"), heavygreen+1pt);
		AddToLegend(format("slope = $%#.2f$", robj.rExec("GetParameter", 1)), heavygreen+1pt);
		AttachLegend(NW, NW);

		NewPad("$\th_y^{*,RN} \ung{\mu rad}$", "$\th_y^{*,RF} \ung{\mu rad}$");
		draw(scale(1e6, 1e6), RootGetObject(f, "optics/p_th_y_RF_vs_th_y_RN"), "d0,eb", red);
		draw(scale(1e6, 1e6), RootGetObject(f, "optics/p_th_y_RF_vs_th_y_RN|pol1"), heavygreen+1pt);
		AddToLegend(format("slope = $%#.2f$", robj.rExec("GetParameter", 1)), heavygreen+1pt);
		AttachLegend(NW, NW);

		NewPad("$\th_y^{*,L} \ung{\mu rad}$", "$\th_y^{*,R} \ung{\mu rad}$");
		draw(scale(1e6, 1e6), RootGetObject(f, "optics/p_th_y_R_vs_th_y_L"), "d0,eb", red);
		draw(scale(1e6, 1e6), RootGetObject(f, "optics/p_th_y_R_vs_th_y_L|pol1"), heavygreen+1pt);
		AddToLegend(format("slope = $%#.2f$", robj.rExec("GetParameter", 1)), heavygreen+1pt);
		AttachLegend(NW, NW);
		
		NewRow();

		NewPad();
		
		NewPad();

		NewPad("$\th_x^{*,L} \ung{\mu rad}$", "$\th_x^{*,R} \ung{\mu rad}$");
		draw(scale(1e6, 1e6), RootGetObject(f, "optics/p_th_x_R_vs_th_x_L"), "d0,eb", blue);
		draw(scale(1e6, 1e6), RootGetObject(f, "optics/p_th_x_R_vs_th_x_L|pol1"), heavygreen+1pt);
		AddToLegend(format("slope = $%#.2f$", robj.rExec("GetParameter", 1)), heavygreen+1pt);
		AttachLegend(NW, NW);

	}
}
