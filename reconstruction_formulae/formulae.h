//----------------------------------------------------------------------------------------------------
// theta_x, single arm, single RP
//----------------------------------------------------------------------------------------------------

Kinematics theta_x_1_arm_one_pot_hit_LF(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_x = - h.x_L_F / env.L_x_L_F;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_x_1_arm_one_pot_hit_LN(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_x = - h.x_L_N / env.L_x_L_N;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_x_1_arm_one_pot_hit_RN(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_x = + h.x_R_N / env.L_x_R_N;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_x_1_arm_one_pot_hit_RF(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_x = + h.x_R_F / env.L_x_R_F;
	return k;
}

//----------------------------------------------------------------------------------------------------
// theta_x, single arm
//----------------------------------------------------------------------------------------------------

Kinematics theta_x_1_arm_hit(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_x_R = + (h.x_R_N / env.L_x_R_N + h.x_R_F / env.L_x_R_F) / 2.;
	k.th_x_L = - (h.x_L_N / env.L_x_L_N + h.x_L_F / env.L_x_L_F) / 2.;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_x_1_arm_angle(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_x_R = + (h.x_R_F - h.x_R_N) / (env.L_x_R_F - env.L_x_R_N);
	k.th_x_L = - (h.x_L_F - h.x_L_N) / (env.L_x_L_F - env.L_x_L_N);
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_x_1_arm_regr(const HitData &h, const Environment &env)
{
	Kinematics k;

	double D_x_L = - env.L_x_L_N * env.v_x_L_F + env.L_x_L_F * env.v_x_L_N;
	k.th_x_L = (env.v_x_L_F * h.x_L_N - env.v_x_L_N * h.x_L_F) / D_x_L;
	//k.vtx_x_L = (env.L_x_L_F * h.x_L_N - env.L_x_L_N * h.x_L_F) / D_x_L;

	double D_x_R = + env.L_x_R_N * env.v_x_R_F - env.L_x_R_F * env.v_x_R_N;
	k.th_x_R = (env.v_x_R_F * h.x_R_N - env.v_x_R_N * h.x_R_F) / D_x_R;
	//k.vtx_x_R = (-env.L_x_R_F * h.x_R_N + env.L_x_R_N * h.x_R_F) / D_x_R;

	return k;
}

//----------------------------------------------------------------------------------------------------
// theta_x, double arm
//----------------------------------------------------------------------------------------------------

Kinematics theta_x_2_arm_LRavg_hit(const HitData &h, const Environment &env)
{
	Kinematics k = theta_x_1_arm_hit(h, env);
	k.th_x = (k.th_x_R + k.th_x_L) / 2.;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_x_2_arm_LRavg_angle(const HitData &h, const Environment &env)
{
	Kinematics k = theta_x_1_arm_angle(h, env);
	k.th_x = (k.th_x_R + k.th_x_L) / 2.;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_x_2_arm_LRavg_regr(const HitData &h, const Environment &env)
{
	Kinematics k = theta_x_1_arm_regr(h, env);
	k.th_x = (k.th_x_R + k.th_x_L) / 2.;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_x_2_arm_regr(const HitData &h, const Environment &env)
{
	Kinematics k;

	double SLL = + env.L_x_L_F*env.L_x_L_F + env.L_x_L_N*env.L_x_L_N + env.L_x_R_N*env.L_x_R_N + env.L_x_R_F*env.L_x_R_F;
	double SLv = - env.L_x_L_F*env.v_x_L_F - env.L_x_L_N*env.v_x_L_N + env.L_x_R_N*env.v_x_R_N + env.L_x_R_F*env.v_x_R_F;
	double Svv = + env.v_x_L_F*env.v_x_L_F + env.v_x_L_N*env.v_x_L_N + env.v_x_R_N*env.v_x_R_N + env.v_x_R_F*env.v_x_R_F;
	double D = (SLL * Svv - SLv * SLv);
	
	double SLx = - env.L_x_L_F*h.x_L_F - env.L_x_L_N*h.x_L_N + env.L_x_R_N*h.x_R_N + env.L_x_R_F*h.x_R_F;
	double Svx = + env.v_x_L_F*h.x_L_F + env.v_x_L_N*h.x_L_N + env.v_x_R_N*h.x_R_N + env.v_x_R_F*h.x_R_F;
		
	k.th_x = (Svv * SLx - SLv * Svx) / D;
	k.vtx_x = (-SLv * SLx + SLL * Svx) / D;

	return k;
}

//----------------------------------------------------------------------------------------------------
// vtx_x, single arm
//----------------------------------------------------------------------------------------------------

Kinematics vtx_x_1_arm_regr(const HitData &h, const Environment &env)
{
	Kinematics k;

	double D_x_L = - env.L_x_L_N * env.v_x_L_F + env.L_x_L_F * env.v_x_L_N;
	k.th_x_L = (env.v_x_L_F * h.x_L_N - env.v_x_L_N * h.x_L_F) / D_x_L;
	k.vtx_x_L = (env.L_x_L_F * h.x_L_N - env.L_x_L_N * h.x_L_F) / D_x_L;

	double D_x_R = + env.L_x_R_N * env.v_x_R_F - env.L_x_R_F * env.v_x_R_N;
	k.th_x_R = (env.v_x_R_F * h.x_R_N - env.v_x_R_N * h.x_R_F) / D_x_R;
	k.vtx_x_R = (-env.L_x_R_F * h.x_R_N + env.L_x_R_N * h.x_R_F) / D_x_R;

	return k;
}
	
//----------------------------------------------------------------------------------------------------
// vtx_x, double arm
//----------------------------------------------------------------------------------------------------

Kinematics vtx_x_2_arm_LRavg_regr(const HitData &h, const Environment &env)
{
	Kinematics k = vtx_x_1_arm_regr(h, env);

	k.vtx_x = (k.vtx_x_L + k.vtx_x_R) / 2.;

	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics vtx_x_2_arm_regr(const HitData &h, const Environment &env)
{
	Kinematics k;

	double SLL = + env.L_x_L_F*env.L_x_L_F + env.L_x_L_N*env.L_x_L_N + env.L_x_R_N*env.L_x_R_N + env.L_x_R_F*env.L_x_R_F;
	double SLv = - env.L_x_L_F*env.v_x_L_F - env.L_x_L_N*env.v_x_L_N + env.L_x_R_N*env.v_x_R_N + env.L_x_R_F*env.v_x_R_F;
	double Svv = + env.v_x_L_F*env.v_x_L_F + env.v_x_L_N*env.v_x_L_N + env.v_x_R_N*env.v_x_R_N + env.v_x_R_F*env.v_x_R_F;
	double D = (SLL * Svv - SLv * SLv);
	
	double SLx = - env.L_x_L_F*h.x_L_F - env.L_x_L_N*h.x_L_N + env.L_x_R_N*h.x_R_N + env.L_x_R_F*h.x_R_F;
	double Svx = + env.v_x_L_F*h.x_L_F + env.v_x_L_N*h.x_L_N + env.v_x_R_N*h.x_R_N + env.v_x_R_F*h.x_R_F;
		
	k.th_x = (Svv * SLx - SLv * Svx) / D;
	k.vtx_x = (-SLv * SLx + SLL * Svx) / D;

	return k;
}

//----------------------------------------------------------------------------------------------------
// theta_x, single arm, single RP
//----------------------------------------------------------------------------------------------------

Kinematics theta_y_1_arm_one_pot_hit_LF(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_y_L = - h.y_L_F / env.L_y_L_F;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_y_1_arm_one_pot_hit_LN(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_y_L = - h.y_L_N / env.L_y_L_N;
	return k;
}

//----------------------------------------------------------------------------------------------------
// theta_y, single arm
//----------------------------------------------------------------------------------------------------

Kinematics theta_y_1_arm_hit(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_y_R = + (h.y_R_N / env.L_y_R_N + h.y_R_F / env.L_y_R_F) / 2.;
	k.th_y_L = - (h.y_L_N / env.L_y_L_N + h.y_L_F / env.L_y_L_F) / 2.;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_y_1_arm_angle(const HitData &h, const Environment &env)
{
	Kinematics k;
	k.th_y_R = + (h.y_R_F - h.y_R_N) / (env.L_y_R_F - env.L_y_R_N);
	k.th_y_L = - (h.y_L_F - h.y_L_N) / (env.L_y_L_F - env.L_y_L_N);
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_y_1_arm_regr(const HitData &h, const Environment &env)
{
	Kinematics k;
	
	//k.th_y_R = + (env.L_y_R_F*h.y_R_F + env.L_y_R_N*h.y_R_N) / (env.L_y_R_F*env.L_y_R_F + env.L_y_R_N*env.L_y_R_N);
	//k.th_y_L = - (env.L_y_L_F*h.y_L_F + env.L_y_L_N*h.y_L_N) / (env.L_y_L_F*env.L_y_L_F + env.L_y_L_N*env.L_y_L_N);

	double D_y_L = - env.L_y_L_N * env.v_y_L_F + env.L_y_L_F * env.v_y_L_N;
	k.th_y_L = (env.v_y_L_F * h.y_L_N - env.v_y_L_N * h.y_L_F) / D_y_L;
	//k.vty_y_L = (env.L_y_L_F * h.y_L_N - env.L_y_L_N * h.y_L_F) / D_y_L;

	double D_y_R = + env.L_y_R_N * env.v_y_R_F - env.L_y_R_F * env.v_y_R_N;
	k.th_y_R = (env.v_y_R_F * h.y_R_N - env.v_y_R_N * h.y_R_F) / D_y_R;
	//k.vty_y_R = (-env.L_y_R_F * h.y_R_N + env.L_y_R_N * h.y_R_F) / D_y_R;


	return k;
}

//----------------------------------------------------------------------------------------------------
// theta_y, double arm
//----------------------------------------------------------------------------------------------------

Kinematics theta_y_2_arm_LRavg_hit(const HitData &h, const Environment &env)
{
	Kinematics k;
	k = theta_y_1_arm_hit(h, env);
	k.th_y = (k.th_y_R + k.th_y_L) / 2.;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_y_2_arm_LRavg_angle(const HitData &h, const Environment &env)
{
	Kinematics k;
	k = theta_y_1_arm_angle(h, env);
	k.th_y = (k.th_y_R + k.th_y_L) / 2.;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_y_2_arm_LRavg_regr(const HitData &h, const Environment &env)
{
	Kinematics k;
	k = theta_y_1_arm_regr(h, env);
	k.th_y = (k.th_y_R + k.th_y_L) / 2.;
	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics theta_y_2_arm_regr(const HitData &h, const Environment &env)
{
	Kinematics k;

	/*
	double SLL = env.L_y_L_F*env.L_y_L_F + env.L_y_L_N*env.L_y_L_N + env.L_y_R_N*env.L_y_R_N + env.L_y_R_F*env.L_y_R_F;
	double SLy =  - env.L_y_L_F*h.y_L_F - env.L_y_L_N*h.y_L_N + env.L_y_R_N*h.y_R_N + env.L_y_R_F*h.y_R_F;

	k.th_y = SLy / SLL;
	*/

	double SLL = + env.L_y_L_F*env.L_y_L_F + env.L_y_L_N*env.L_y_L_N + env.L_y_R_N*env.L_y_R_N + env.L_y_R_F*env.L_y_R_F;
	double SLv = - env.L_y_L_F*env.v_y_L_F - env.L_y_L_N*env.v_y_L_N + env.L_y_R_N*env.v_y_R_N + env.L_y_R_F*env.v_y_R_F;
	double Svv = + env.v_y_L_F*env.v_y_L_F + env.v_y_L_N*env.v_y_L_N + env.v_y_R_N*env.v_y_R_N + env.v_y_R_F*env.v_y_R_F;
	double D = (SLL * Svv - SLv * SLv);
	
	double SLy = - env.L_y_L_F*h.y_L_F - env.L_y_L_N*h.y_L_N + env.L_y_R_N*h.y_R_N + env.L_y_R_F*h.y_R_F;
	double Svy = + env.v_y_L_F*h.y_L_F + env.v_y_L_N*h.y_L_N + env.v_y_R_N*h.y_R_N + env.v_y_R_F*h.y_R_F;
		
	k.th_y = (Svv * SLy - SLv * Svy) / D;
	k.vtx_y = (-SLv * SLy + SLL * Svy) / D;

	return k;
}

//----------------------------------------------------------------------------------------------------
// vtx_y, single arm
//----------------------------------------------------------------------------------------------------

Kinematics vtx_y_1_arm_regr(const HitData &h, const Environment &env)
{
	Kinematics k;

	double D_y_L = - env.L_y_L_N * env.v_y_L_F + env.L_y_L_F * env.v_y_L_N;
	k.th_y_L = (env.v_y_L_F * h.y_L_N - env.v_y_L_N * h.y_L_F) / D_y_L;
	k.vtx_y_L = (env.L_y_L_F * h.y_L_N - env.L_y_L_N * h.y_L_F) / D_y_L;

	double D_y_R = + env.L_y_R_N * env.v_y_R_F - env.L_y_R_F * env.v_y_R_N;
	k.th_y_R = (env.v_y_R_F * h.y_R_N - env.v_y_R_N * h.y_R_F) / D_y_R;
	k.vtx_y_R = (-env.L_y_R_F * h.y_R_N + env.L_y_R_N * h.y_R_F) / D_y_R;

	return k;
}
	
//----------------------------------------------------------------------------------------------------
// vtx_y, double arm
//----------------------------------------------------------------------------------------------------

Kinematics vtx_y_2_arm_LRavg_regr(const HitData &h, const Environment &env)
{
	Kinematics k = vtx_y_1_arm_regr(h, env);

	k.vtx_y = (k.vtx_y_L + k.vtx_y_R) / 2.;

	return k;
}

//----------------------------------------------------------------------------------------------------

Kinematics vtx_y_2_arm_regr(const HitData &h, const Environment &env)
{
	Kinematics k;

	double SLL = + env.L_y_L_F*env.L_y_L_F + env.L_y_L_N*env.L_y_L_N + env.L_y_R_N*env.L_y_R_N + env.L_y_R_F*env.L_y_R_F;
	double SLv = - env.L_y_L_F*env.v_y_L_F - env.L_y_L_N*env.v_y_L_N + env.L_y_R_N*env.v_y_R_N + env.L_y_R_F*env.v_y_R_F;
	double Svv = + env.v_y_L_F*env.v_y_L_F + env.v_y_L_N*env.v_y_L_N + env.v_y_R_N*env.v_y_R_N + env.v_y_R_F*env.v_y_R_F;
	double D = (SLL * Svv - SLv * SLv);
	
	double SLx = - env.L_y_L_F*h.y_L_F - env.L_y_L_N*h.y_L_N + env.L_y_R_N*h.y_R_N + env.L_y_R_F*h.y_R_F;
	double Svx = + env.v_y_L_F*h.y_L_F + env.v_y_L_N*h.y_L_N + env.v_y_R_N*h.y_R_N + env.v_y_R_F*h.y_R_F;
		
	k.th_y = (Svv * SLx - SLv * Svx) / D;
	k.vtx_y = (-SLv * SLx + SLL * Svx) / D;

	return k;
}
