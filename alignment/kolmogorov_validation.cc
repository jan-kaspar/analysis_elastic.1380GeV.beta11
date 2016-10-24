#include "kolmogorov.h"

#include "TFile.h"

int main()
{
	TFile *f_out = TFile::Open("kolmogorov_validation.root", "recreate");

	CumulativeHistogram h1, h2;

	h1.Fill(1.);
	h1.Fill(2.);
	h1.Fill(5.);
	h1.Fill(2.);
	h1.GetGraph(false)->Write("h1");

	h2.Fill(1.5);
	h2.Fill(2.);
	h2.Fill(4.);
	h2.Fill(4.5);
	h2.GetGraph(false)->Write("h2");

	TGraph *g_diff = new TGraph();

	printf("test: %f\n", KolmogorovTest(h1, h2, false, g_diff));

	g_diff->Write("g_diff");

	delete f_out;
	return 0;
}
