void Analysis::BuildCuts()
{
	N_cuts = 8;

	// cut structure:
	//	| a*qa + b*qb + c| < n_si * si

	// a: th_x_R, b: th_x_L
	cqaN[1] = "#theta_{x}^{R}"; cqbN[1] = "#theta_{x}^{L}";
	cca[1] = -cut1_a;
	ccb[1] = 1.;
	ccc[1] = cut1_c;
	csi[1] = cut1_si;
	cuts.push_back(1);
	
	// a: th_y_R, b: th_y_L
	cqaN[2] = "#theta_{y}^{R}"; cqbN[2] = "#theta_{y}^{L}";
	cca[2] = -cut2_a;
	ccb[2] = 1.;
	ccc[2] = cut2_c;
	csi[2] = cut2_si;
	cuts.push_back(2);

	// a: th_x_R, b: vtx_x_R
	cqaN[3] = "#theta_{x}^{R}"; cqbN[3] = "vtx_{x}^{R}";
	cca[3] = -cut3_a;
	ccb[3] = 1.;
	ccc[3] = cut3_b;
	csi[3] = cut3_si;
	//cuts.push_back(3);
	
	// a: th_x_L, b: vtx_x_L
	cqaN[4] = "#theta_{x}^{L}"; cqbN[4] = "vtx_{x}^{L}";
	cca[4] = -cut4_a;
	ccb[4] = 1.;
	ccc[4] = cut4_b;
	csi[4] = cut4_si;
	//cuts.push_back(4);
	
	// a: y_R_N, b: y_R_F - y_R_N
	cqaN[5] = "y^{R,N}"; cqbN[5] = "y^{R,F} - y^{R,N}";
	cca[5] = -cut5_a;
	ccb[5] = 1.;
	ccc[5] = cut5_b;
	csi[5] = cut5_si;
	// TODO: uncomment
	//cuts.push_back(5);
	
	// a: y_L_N, b: y_L_F - y_L_N
	cqaN[6] = "y^{L,N}"; cqbN[6] = "y^{L,F} - y^{L,N}";
	cca[6] = -cut6_a;
	ccb[6] = 1.;
	ccc[6] = cut6_b;
	csi[6] = cut6_si;
	// TODO: uncomment
	//cuts.push_back(6);
	
	// a: th_x, b: vtx_x_R - vtx_x_L
	cqaN[7] = "#theta_{x}"; cqbN[7] = "vtx_{x}^{R} - vtx_{x}^{L}";
	cca[7] = -cut7_a;
	ccb[7] = 1.;
	ccc[7] = cut7_c;
	csi[7] = cut7_si;
	cuts.push_back(7);
	
	// a: th_y, b: vtx_y_R - vtx_y_L
	cqaN[8] = "#theta_{y}"; cqbN[8] = "vtx_{y}^{R} - vtx_{y}^{L}";
	cca[8] = -cut8_a;
	ccb[8] = 1.;
	ccc[8] = cut8_c;
	csi[8] = cut8_si;
}

//----------------------------------------------------------------------------------------------------

bool Analysis::EvaluateCuts(const HitData &h, const Kinematics &k, CutData &cd) const
{
	cd.cqa[1] = k.th_x_R;	cd.cqb[1] = k.th_x_L;
	cd.cqa[2] = k.th_y_R;	cd.cqb[2] = k.th_y_L;
	cd.cqa[3] = k.th_x_R;	cd.cqb[3] = k.vtx_x_R;
	cd.cqa[4] = k.th_x_L;	cd.cqb[4] = k.vtx_x_L;
	cd.cqa[5] = h.y_R_N;	cd.cqb[5] = h.y_R_F - h.y_R_N;
	cd.cqa[6] = h.y_L_N;	cd.cqb[6] = h.y_L_F - h.y_L_N;
	cd.cqa[7] = k.th_x;		cd.cqb[7] = k.vtx_x_R - k.vtx_x_L;
	cd.cqa[8] = k.th_y;		cd.cqb[8] = k.vtx_y_R - k.vtx_y_L;

	for (unsigned int ci = 1; ci <= N_cuts; ++ci)
	{
		cd.cv[ci] = cca[ci]*cd.cqa[ci] + ccb[ci]*cd.cqb[ci] + ccc[ci];
		cd.ct[ci] = (fabs(cd.cv[ci]) <= n_si * csi[ci]);
		//printf("cut %u: |%+E| < %E * %E <==> %i\n", ci, cd.cv[ci], n_si, csi[ci], cd.ct[ci]);
	}

	bool select = true;	// and between all cuts
	for (unsigned int ci = 0; ci < cuts.size(); ci++)
	{
		select &= cd.ct[cuts[ci]];
	}

	return select;
}
