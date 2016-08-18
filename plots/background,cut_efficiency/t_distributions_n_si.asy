import root;
import pad_layout;

string datasets[] = { "DS2", "DS4" };

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b -- 56t", "45t -- 56b" };

string topDir = "../../";

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void PlotRatio(rObject o, rObject r, pen p)
{
	int N = r.iExec("GetNbinsX");
	for (int i = 1; i <= N; ++i)
	{
		real c = r.rExec("GetBinCenter", i);
		real v = o.rExec("GetBinContent", i);
		real v_r = r.rExec("GetBinContent", i);

		real ratio = (v_r != 0) ? (v / v_r - 1.) * 100.: 0.;
		real ratio_unc = 0. * 100;

		draw(Scale((c, ratio)), mCi+p);
		draw(Scale((c, ratio-ratio_unc))--Scale((c, ratio+ratio_unc)), p);
	}
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

string combinations[];
pen comb_pens[];

combinations.push("n_si=2"); comb_pens.push(black);
combinations.push("n_si=3"); comb_pens.push(red);
combinations.push("n_si=4"); comb_pens.push(blue);
combinations.push("n_si=5"); comb_pens.push(heavygreen);
combinations.push("n_si=6"); comb_pens.push(magenta);


string ref_comb = "n_si=4";

xSizeDef = 10cm;
xTicksDef=LeftTicks(0.2, 0.1);

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	for (int dgi : diagonals.keys)
	{
		NewRow();

		string dataset = datasets[dsi];
		string diagonal = diagonals[dgi];

		NewPad(false);
		label(replace("{\SetFontSizesXX\vbox{\hbox{"+dataset+"}\hbox{"+dgn_labels[dgi]+"}}}", "_", "\_"));

		NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t\ung{mb/GeV^2}$");
		scale(Linear, Log);
		for (int ci : combinations.keys)
		{
			string label = replace(combinations[ci], "_", "\_");
			if (combinations[ci] == ref_comb)
				label += " (reference)";
		
			string f = topDir+dataset+"/background_study/"+combinations[ci]+"/distributions_"+diagonal+".root";
			draw(rGetObj(f, "normalization/eb/h_t_normalized"), "eb",
				comb_pens[ci], label);
		}
		
		frame fLegend = BuildLegend();
		
		limits((0, 1e-3), (1.4, 1e3), Crop);
		
		NewPad(false);
		attach(fLegend);
		
		//--------------------
		
		string ref_f = topDir+dataset+"/background_study/"+ref_comb+"/distributions_"+diagonal+".root";
		rObject ref_o = rGetObj(ref_f, "normalization/eb/h_t_normalized");
		
		//--------------------

		NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t$: (test - ref) / ref$\ung{\%}$");
		for (int ci : combinations.keys)
		{
			//if (combinations[ci] == ref_comb)
			//	continue;
		
			string f = topDir+dataset+"/background_study/"+combinations[ci]+"/distributions_"+diagonal+".root";
			rObject o = rGetObj(f, "normalization/eb/h_t_normalized");
			PlotRatio(o, ref_o, comb_pens[ci]);
		}
		
		limits((0, -5), (1.4, 10), Crop);
		
		//--------------------
		
		/*
		NewPad("$|t|\ung{GeV^2}$", "$\d\sigma/\d t$: (test - ref) / ref$\ung{\%}$", xTicks=LeftTicks());
		for (int ci : combinations.keys)
		{
			//if (combinations[ci] == ref_comb)
			//	continue;
		
			string f = topDir+dataset+"/background_study/"+combinations[ci]+"/distributions_"+diagonal+".root";
			rObject o = rGetObj(f, "normalization/eb/h_t_normalized");
			PlotRatio(o, ref_o, comb_pens[ci]);
		}
		
		limits((0, -5), (0.3, 5), Crop);
		*/
	}
}
