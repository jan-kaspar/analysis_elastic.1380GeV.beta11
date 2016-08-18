void test()
{
	TChain *ch = new TChain("ntuple");
	ch->Add("ntuple_sim_hits.root");

	//ch->Draw("y:z", "abs(z - -220.000E3) < 100");
	//ch->Draw("y:z", "abs(z - -214.628E3) < 100");
	//ch->Draw("y:z", "abs(z - 214.628E3) < 100");
	ch->Draw("y:z", "abs(z - 220.000E3) < 100");
	
	//ch->Draw("y:z", "(z > 219.95E3) && (z < 220.05E3) && (Event == 616)");
}
