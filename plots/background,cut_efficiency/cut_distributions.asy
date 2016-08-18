import root;
import pad_layout;

string datasets[] = { "DS1" };

string diagonals[] = { "45b_56t", "45t_56b" };
string dgn_labels[] = { "45b -- 56t", "45t -- 56b" };

string topDir = "../../";

xSizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

string dataset;
string diagonal;

void MakeComparison(string quantity, real xscale, string unit, string obj, real xlimit, real sigma,
	real xStep, real xstep,
	string combinations[], pen comb_pens[])
{
	NewPad(quantity+"$\ung{"+unit+"}$", "", xTicks=LeftTicks(xStep, xstep));
	scale(Linear, Log(true));
	for (int ci : combinations.keys)
	{
		string f = topDir+dataset+"/background_study/"+combinations[ci]+"/distributions_"+diagonal+".root";
		draw(scale(xscale, 1), rGetObj(f, "elastic cuts/"+obj), "vl", 
			comb_pens[ci], replace(combinations[ci], "_", "\_"));	
	}

	yaxis(XEquals(-4*sigma, false), dashed);
	yaxis(XEquals(+4*sigma, false), dashed);

	xlimits(-xlimit, +xlimit, Crop);
	AttachLegend(quantity, NW, NE);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
{
	for (int dgi : diagonals.keys)
	{
		NewRow();

		dataset = datasets[dsi];
		diagonal = diagonals[dgi];

		NewPad(false);
		label(replace("{\SetFontSizesXX\vbox{\hbox{"+dataset+"}\hbox{"+dgn_labels[dgi]+"}}}", "_", "\_"));

		string combinations[];
		pen comb_pens[];
		
		combinations.push("no_cuts"); comb_pens.push(gray);
		combinations.push("cuts:2"); comb_pens.push(red);
		combinations.push("cuts:2,5"); comb_pens.push(blue);
		combinations.push("cuts:2,5,6"); comb_pens.push(magenta);
		combinations.push("cuts:2,5,6,7"); comb_pens.push(heavygreen);
		//combinations.push("cuts:2,5,6,7,3"); comb_pens.push(cyan);
		//combinations.push("cuts:2,5,6,7,3,4"); comb_pens.push(orange);
		
		MakeComparison("discriminator 1: $\De\th_x^*$", 1e6, "\mu rad", "cut 1/h_cq1", 400, 9.4, 100, 20, combinations, comb_pens);
		
		//--------------------
		//NewPage();
		
		string combinations[];
		pen comb_pens[];
		
		combinations.push("no_cuts"); comb_pens.push(gray);
		combinations.push("cuts:1"); comb_pens.push(red);
		combinations.push("cuts:1,5"); comb_pens.push(blue);
		combinations.push("cuts:1,5,6"); comb_pens.push(magenta);
		combinations.push("cuts:1,5,6,7"); comb_pens.push(heavygreen);
		//combinations.push("cuts:1,5,6,7,3"); comb_pens.push(cyan);
		//combinations.push("cuts:1,5,6,7,3,4"); comb_pens.push(orange);
		
		MakeComparison("discriminator 2: $\De\th_y^*$", 1e6, "\mu rad", "cut 2/h_cq2", 200, 3.3, 50, 10, combinations, comb_pens);
		
		//--------------------
		//NewPage();
		
		string combinations[];
		pen comb_pens[];
		
		combinations.push("no_cuts"); comb_pens.push(gray);
		combinations.push("cuts:1"); comb_pens.push(red);
		combinations.push("cuts:1,2"); comb_pens.push(blue);
		combinations.push("cuts:1,2,6"); comb_pens.push(magenta);
		combinations.push("cuts:1,2,6,7"); comb_pens.push(heavygreen);
		
		MakeComparison("discriminator 5: $\De^{F-N} y^{R}\hbox{ vs. }y^{RN}$", 1e3, "\mu m", "cut 5/h_cq5", 500, 18., 100, 20, combinations, comb_pens);
		
		//--------------------
		//NewPage();
		
		string combinations[];
		pen comb_pens[];
		
		combinations.push("no_cuts"); comb_pens.push(gray);
		combinations.push("cuts:1"); comb_pens.push(red);
		combinations.push("cuts:1,2"); comb_pens.push(blue);
		combinations.push("cuts:1,2,5"); comb_pens.push(magenta);
		combinations.push("cuts:1,2,5,7"); comb_pens.push(heavygreen);
		
		MakeComparison("discriminator 6: $\De^{F-N} y^{L}\hbox{ vs. }y^{LN}$", 1e3, "\mu m", "cut 6/h_cq6", 500, 18., 100, 20, combinations, comb_pens);
		
		//--------------------
		//NewPage();
		
		string combinations[];
		pen comb_pens[];
		
		combinations.push("no_cuts"); comb_pens.push(gray);
		combinations.push("cuts:1"); comb_pens.push(red);
		combinations.push("cuts:1,2"); comb_pens.push(blue);
		combinations.push("cuts:1,2,5"); comb_pens.push(magenta);
		combinations.push("cuts:1,2,5,6"); comb_pens.push(heavygreen);
		//combinations.push("cuts:1,2,5,6,3"); comb_pens.push(cyan);
		//combinations.push("cuts:1,2,5,6,3,4"); comb_pens.push(orange);
		
		MakeComparison("discriminator 7: $\De x^*$", 1e3, "\mu m", "cut 7/h_cq7", 500, 8.7, 100, 20, combinations, comb_pens);
	}
}
