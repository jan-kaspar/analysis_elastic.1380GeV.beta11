#include <vector>
#include <cmath>
#include <cstdio>

using namespace std;

int main()
{
	double v_x_L_N = -3.211961106, v_x_R_N = -3.260793214;		// 1
	double v_x_L_F = -2.943772831, v_x_R_F = -3.003100052;		// 1

	double L_x_L_N = 6.907815061, L_x_R_N = 7.123123432;		// m
	double L_x_L_F = 4.658536578, L_x_R_F = 4.912747084;		// m

	double v_y_L_N = -3.241671983, v_y_R_N = -3.188960453;		// 1
	double v_y_L_F = -3.51398971, v_y_R_F = -3.45003027;		// 1

	double L_y_L_N = 20.041446776, L_y_R_N = 19.452649165;		// m
	double L_y_L_F = 20.067865623, L_y_R_F = 19.360612796;		// m

	vector< vector<double> > modes;

	double rel_unc = 0.07;

	// v_x_L_N, L_x_L_N, v_y_L_N, L_y_L_N, v_x_L_F, L_x_L_F, v_y_L_F, L_y_L_F
	vector<double> x_mode;
	x_mode.push_back(rel_unc * fabs(v_x_L_N));
	x_mode.push_back(rel_unc * fabs(L_x_L_N));
	x_mode.push_back(0.);
	x_mode.push_back(0.);
	x_mode.push_back(rel_unc * fabs(v_x_L_F));
	x_mode.push_back(rel_unc * fabs(L_x_L_F));
	x_mode.push_back(0.);
	x_mode.push_back(0.);
	modes.push_back(x_mode);
	
	vector<double> y_mode;
	y_mode.push_back(0.);
	y_mode.push_back(0.);
	y_mode.push_back(rel_unc * fabs(v_y_L_N));
	y_mode.push_back(rel_unc * fabs(L_y_L_N));
	y_mode.push_back(0.);
	y_mode.push_back(0.);
	y_mode.push_back(rel_unc * fabs(v_y_L_F));
	y_mode.push_back(rel_unc * fabs(L_y_L_F));
	modes.push_back(y_mode);

	// initiate matrix
	unsigned int N = 8;
	vector< vector<double> > matrix;
	for (unsigned int i = 0; i < N; i++)
	{
		vector<double> v;
		for (unsigned int j = 0; j < N; j++)
			v.push_back(0.);

		matrix.push_back(v);
	}

	// build matrix
	for (unsigned int mi = 0; mi < modes.size(); mi++)
	{
		for (unsigned int i = 0; i < N; i++)
			for (unsigned int j = 0; j < N; j++)
			{
				matrix[i][j] += modes[mi][i] * modes[mi][j];
			}
	}

	// print matrix
	for (unsigned int i = 0; i < N; i++)
	{
		for (unsigned int j = 0; j < N; j++)
		{
			printf("%.5E,\t", matrix[i][j]);
		}
		printf("\n");
	}
}
