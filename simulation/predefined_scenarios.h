int SetPredefinedScenario(double /*th_y_sign*/, Environment &/*env_nom*/, Environment &env_sim,
		Analysis &/*anal_id*/, Analysis &/*anal_rec*/, const string &scenario)
{
	// uncertainties
	double si_de_p = env_sim.si_de_p;
	
	double si_v_x_N = sqrt(env_sim.opt_cov(0, 0));
	double si_v_x_F = sqrt(env_sim.opt_cov(4, 4));
                                                   
	double si_L_x_N = sqrt(env_sim.opt_cov(1, 1)) * 1E3;	// m --> mm
	double si_L_x_F = sqrt(env_sim.opt_cov(5, 5)) * 1E3;
                                                   
	double si_L_y_N = sqrt(env_sim.opt_cov(3, 3)) * 1E3;
	double si_L_y_F = sqrt(env_sim.opt_cov(7, 7)) * 1E3;

	double si_de_x = env_sim.si_de_x;
	double si_de_y = env_sim.si_de_y;
	double si_tilt = env_sim.si_tilt;

	// ideal situation
	if (scenario.compare("none") == 0)
	{
		return 0;
	}

#if 0
	// various values of si_th_x and si_th_y (to study unsmearing correction)
	// TODO: values of si_th_x and si_th_y
	if (scenario.compare("nom.si_th_x=0.47") == 0)
	{
		env_nom.si_th_x = env_sim.si_th_x = 0.47E-6;
		return 10;
	}
	
	if (scenario.compare("nom.si_th_y=0.65") == 0)
	{
		env_nom.si_th_y = env_sim.si_th_y = anal_id.si_th_y_1arm = anal_rec.si_th_y_1arm = 0.65E-6;
		return 10;
	}

	// error of si_th_x,y (for acceptance correction)
	if (scenario.compare("si_th_y-1_th_x-1") == 0)
	{
		anal_rec.si_th_y_1arm +=	-1. * anal_rec.si_th_y_1arm_unc;
		anal_rec.si_th_x_1arm +=	-1. * anal_rec.si_th_x_1arm_unc;
		return 0;
	}
#endif
	
	
	// vertical misalignment
	if (scenario.compare("de_y+1+1") == 0)
	{
		env_sim.de_y_R_N += +1.*si_de_y;
		env_sim.de_y_R_F += +1.*si_de_y;
		return 0;
	}

	if (scenario.compare("de_y+1-1") == 0)
	{
		env_sim.de_y_R_N += +1.*si_de_y;
		env_sim.de_y_R_F += -1.*si_de_y;
		return 0;
	}
	
	if (scenario.compare("de_y-1-1") == 0)
	{
		env_sim.de_y_R_N += -1.*si_de_y;
		env_sim.de_y_R_F += -1.*si_de_y;
		return 0;
	}

	// horizontal misalignment
	if (scenario.compare("de_x+1+1") == 0)
	{
		env_sim.de_x_R_N += +1.*si_de_x;
		env_sim.de_x_R_F += +1.*si_de_x;
		return 0;
	}
	
	if (scenario.compare("de_x-1+1") == 0)
	{
		env_sim.de_x_R_N += -1.*si_de_x;
		env_sim.de_x_R_F += +1.*si_de_x;
		return 0;
	}
	
	if (scenario.compare("de_x+1-1") == 0)
	{
		env_sim.de_x_R_N += +1.*si_de_x;
		env_sim.de_x_R_F += -1.*si_de_x;
		return 0;
	}
	
	// tilt misalignment
	if (scenario.compare("tilt+1+1") == 0)
	{
		env_sim.tilt_R_N += +1.*si_tilt;
		env_sim.tilt_R_F += +1.*si_tilt;
		return 0;
	}
	
	if (scenario.compare("tilt+1-1") == 0)
	{
		env_sim.tilt_R_N += +1.*si_tilt;
		env_sim.tilt_R_F += -1.*si_tilt;
		return 0;
	}
	
	// optics vertical error
	if (scenario.compare("L_y+1+1") == 0)
	{
		env_sim.L_y_R_N += +1.*si_L_y_N;
		env_sim.L_y_R_F += +1.*si_L_y_F;
		return 0;
	}
	
	if (scenario.compare("L_y-1+1") == 0)
	{
		env_sim.L_y_R_N += -1.*si_L_y_N;
		env_sim.L_y_R_F += +1.*si_L_y_F;
		return 0;
	}
	
	// optics horizontal error
	if (scenario.compare("L_x+1+1") == 0)
	{
		env_sim.L_x_R_N += +1.*si_L_x_N;
		env_sim.L_x_R_F += +1.*si_L_x_F;
		return 0;
	}
	
	if (scenario.compare("L_x+1-1") == 0)
	{
		env_sim.L_x_R_N += +1.*si_L_x_N;
		env_sim.L_x_R_F += -1.*si_L_x_F;
		return 0;
	}
	if (scenario.compare("v_x+1+1") == 0)
	{
		env_sim.v_x_R_N += +1.*si_v_x_N;
		env_sim.v_x_R_F += +1.*si_v_x_F;
		return 0;
	}
	
	if (scenario.compare("v_x+1-1") == 0)
	{
		env_sim.v_x_R_N += +1.*si_v_x_N;
		env_sim.v_x_R_F += -1.*si_v_x_F;
		return 0;
	}
	
	// optics modes
	if (scenario.compare("opt+1m1") == 0)
	{
		env_sim.ApplyOpticsPerturbationMode(0, 1.);
		return 0;
	}
	
	if (scenario.compare("opt+1m2") == 0)
	{
		env_sim.ApplyOpticsPerturbationMode(1, 1.);
		return 0;
	}
	
	if (scenario.compare("opt+1m3") == 0)
	{
		env_sim.ApplyOpticsPerturbationMode(2, 1.);
		return 0;
	}
	
	if (scenario.compare("opt+1m4") == 0)
	{
		env_sim.ApplyOpticsPerturbationMode(3, 1.);
		return 0;
	}
	
	if (scenario.compare("opt+1m5") == 0)
	{
		env_sim.ApplyOpticsPerturbationMode(4, 1.);
		return 0;
	}
	
	if (scenario.compare("opt+1m6") == 0)
	{
		env_sim.ApplyOpticsPerturbationMode(5, 1.);
		return 0;
	}
	
	if (scenario.compare("opt+1m7") == 0)
	{
		env_sim.ApplyOpticsPerturbationMode(6, 1.);
		return 0;
	}
	
	if (scenario.compare("opt+1m8") == 0)
	{
		env_sim.ApplyOpticsPerturbationMode(7, 1.);
		return 0;
	}
	
	// beam momentum error
	if (scenario.compare("de_p+1") == 0)
	{
		env_sim.p += +1.*si_de_p;
		return 0;
	}
	
/*
	// error of normalization
	if (scenario.compare("norm+1") == 0)
	{
		anal_rec.full_norm_corr += +1. * anal_rec.full_norm_corr_unc;
		return 0;
	}

	if (scenario.compare("norm-1") == 0)
	{
		anal_rec.full_norm_corr += -1. * anal_rec.full_norm_corr_unc;
		return 0;
	}
*/

	return 1; 
}
