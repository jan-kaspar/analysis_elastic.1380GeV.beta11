import root;
import pad_layout;

string topDir = "../../";

string datasets[] = {
	"DS1",
//	"DS3"
};

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b -- 56t", "45t -- 56b" };
pen dgn_pens[] = { blue, red };

string arms[] = { "_L", "_R", "" };
string arm_ss[] = { "L", "R", "" };
string arm_labels[] = { "left", "right", "double" };

xSizeDef = 12cm;
//xTicksDef = LeftTicks(Step=1, step=0.5);

//TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

NewPad(false);
for (int ai : arms.keys)
{
	NewPad(false);
	label("{\SetFontSizesXX " + arm_labels[ai] + " arm}");
}

//----------------------------------------------------------------------------------------------------

frame f_legend;

for (int dsi : datasets.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX " + datasets[dsi] + "}");

	for (int ai : arms.keys)
	{
		NewPad("$\th_y^{*"+arm_ss[ai]+"}\ung{\mu rad}$", "mean of~$\th_x^{*"+arm_ss[ai]+"}$");
		for (int dgi : diagonals.keys)
		{
			TF1_x_min = -inf;
			TF1_x_max = +inf;

			string f = topDir+datasets[dsi]+"/distributions_"+diagonals[dgi]+".root";
			
			draw(scale(1e6, 1e6), RootGetObject(f, "selected - angles/p_th_x"+arms[ai]+"_vs_th_y"+arms[ai]),
				"eb,d0", StdPen(dgi+1), dgn_labels[dgi]);
			
			RootObject fit = RootGetObject(f, "selected - angles/p_th_x"+arms[ai]+"_vs_th_y"+arms[ai]+"|pol1", error=false);
			if (!fit.valid)
				continue;

			draw(scale(1e6, 1e6), fit, black+dashed, (dgi == 0) ? "single-RP fit" : "");
	
			real a = fit.rExec("GetParameter", 1);
			real a_u = fit.rExec("GetParError", 1);
	
			real x = (dgi == 0) ? +300 : -300;
			label(format("slope = $%.3f$", a) + format("$\pm %.3f$", a_u), (x, 70));
		}

		/*
		string ff = topDir + "overall_alignment/fit_mean_th_x_vs_th_y.root";
		//draw(scale(1e6, 1e6), RootGetObject(ff, datasets[dsi] + "/" + arm_labels[ai] + "/g_mean_th_x_vs_th_y"), "p", mCi+3pt);		
		TF1_x_min = -90e-6;
		TF1_x_max = +90e-6;

		RootObject fit = RootGetObject(ff, datasets[dsi] + "/" + arm_labels[ai] + "/pol1");
		draw(scale(1e6, 1e6), fit, heavygreen+2pt, "double-RP fit");

		real a = fit.rExec("GetParameter", 1);
		real a_u = fit.rExec("GetParError", 1);
		real chisq = fit.rExec("GetChisquare");
		int ndf = fit.iExec("GetNDF");

		real x = 0;
		label(format("$\ch^2 / \hbox{ndf} = %.3f$", chisq) + format("$/ %i$", ndf) + format("$ = %.3f$", chisq/ndf), (x, +2.5), heavygreen);
		label(format("slope = $%.3f$", a) + format("$\pm %.3f$", a_u), (x, +2), heavygreen);
		*/

		limits((-500, -100), (+500, +100), Crop);

		f_legend = BuildLegend();
	}

	NewPad(false);
	add(f_legend);
}

GShipout(vSkip=0mm);
