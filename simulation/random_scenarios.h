int GenerateRandomClass(double /*th_y_sign*/, Environment &env_sim, Analysis &/*anal_id*/, Analysis &anal_rec, const string &c)
{
	// experimental angular resolutions
	if (c.compare("res") == 0)
	{
		anal_rec.si_th_y_1arm +=	gRandom->Gaus() * anal_rec.si_th_y_1arm_unc;
		anal_rec.si_th_x_1arm_L +=	gRandom->Gaus() * anal_rec.si_th_x_1arm_unc;
		anal_rec.si_th_x_1arm_R +=	gRandom->Gaus() * anal_rec.si_th_x_1arm_unc;
		return 0;
	}

	// beam energies
	if (c.compare("bmom") == 0)
	{
		double si_de_p = env_sim.si_de_p;
		env_sim.p_L += gRandom->Gaus() * si_de_p;
		env_sim.p_R += gRandom->Gaus() * si_de_p;
		env_sim.p = (env_sim.p_L + env_sim.p_R) / 2.;
		return 0;
	}

	// misalignments
	if (c.compare("misal") == 0)
	{
		double si_de_x = env_sim.si_de_x;
		double si_de_y = env_sim.si_de_y;
		double si_tilt = env_sim.si_tilt;

		env_sim.de_x_L_F += gRandom->Gaus() * si_de_x; env_sim.de_y_L_F += gRandom->Gaus() * si_de_y; env_sim.tilt_L_F += gRandom->Gaus() * si_tilt;
		env_sim.de_x_L_N += gRandom->Gaus() * si_de_x; env_sim.de_y_L_N += gRandom->Gaus() * si_de_y; env_sim.tilt_L_N += gRandom->Gaus() * si_tilt;
		env_sim.de_x_R_N += gRandom->Gaus() * si_de_x; env_sim.de_y_R_N += gRandom->Gaus() * si_de_y; env_sim.tilt_R_N += gRandom->Gaus() * si_tilt;
		env_sim.de_x_R_F += gRandom->Gaus() * si_de_x; env_sim.de_y_R_F += gRandom->Gaus() * si_de_y; env_sim.tilt_R_F += gRandom->Gaus() * si_tilt;

		return 0;
	}
	
	// optics
	if (c.compare("optics") == 0)
	{
		TVectorD de(16);
		env_sim.ApplyRandomOpticsPerturbations(de);

		return 0;
	}

	// normalization
/*
	if (c.compare("norm") == 0)
	{
		anal_rec.full_norm_corr += gRandom->Gaus() * anal_rec.full_norm_corr_unc;
		return 0;
	}
*/

	return 1;
}

//----------------------------------------------------------------------------------------------------

int GenerateRandomScenario(double th_y_sign, Environment &env_sim, Analysis &anal_id, Analysis &anal_rec, const string &list)
{
	string::size_type start = 0, end = 0;
	while (end != string::npos)
	{
		end = list.find(",", start);
		string c = (end == string::npos) ? list.substr(start) : list.substr(start, end - start);
		int result = GenerateRandomClass(th_y_sign, env_sim, anal_id, anal_rec, c);
		if (result)
		{
			printf("ERROR: unknown random class `%s'\n", c.c_str());
			return 1;
		}
		start = end + 1;
	}

	return 0; 
}
