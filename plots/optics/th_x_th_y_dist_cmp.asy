import root;
import pad_layout;

string topDir = "../../";

string datasets[] = { "DS1", "DS3" };
string dataset_labels[] = { "DS2, matched opt.", "DS3, matched opt." };

//string arms[] = { "45_top", "56_top", "45_bottom", "56_bottom" };
//string arm_labels[] = { "45 top", "56 top", "45 bottom", "56 bottom" };


xSizeDef = 8cm;
//xTicksDef = LeftTicks(0.02, 0.01);


for (int dsi : datasets.keys)
{
	NewRow();

	NewPad(false);
	label("{\SetFontSizesXX\vbox{\hbox{" + dataset_labels[dsi] + "}} }");

	string f = topDir+datasets[dsi]+"/scale_th_x.root";

	// --------------------
	NewPad("$\th_{x,y}^*\ung{\mu rad}$");

	draw(rGetObj(f, "global/fit/h_th_x"), "eb,d0", red, "$\th_x^*$");
	draw(rGetObj(f, "global/fit/h_th_y"), "eb,d0", blue, "$\th_y^*$");

	AttachLegend();
}
