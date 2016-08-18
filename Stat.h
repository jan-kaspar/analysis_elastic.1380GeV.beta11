#include <vector>

using namespace std;


struct Stat
{
	unsigned int dim;
	double S1;
	vector<double> Sv;
	vector< vector<double> > Svv;

	Stat(unsigned int _dim = 1)
	{
		dim = _dim;

		S1 = 0.;
		for (unsigned int i = 0; i < dim; i++) {
			Sv.push_back(0);

			vector<double> temp;
			for (unsigned int j = 0; j < dim; j++) {
				temp.push_back(0);
			}
			Svv.push_back(temp);
		}
	}

	void Fill(const vector<double> &v)
	{
		S1 += 1.;
		for (unsigned int i = 0; i < dim; i++) {
			Sv[i] += v[i];

			for (unsigned int j = 0; j < dim; j++) {
				Svv[i][j] += v[i] * v[j];
			}
		}
	}

	void Fill(const TVectorD &v)
	{
		S1 += 1.;
		for (unsigned int i = 0; i < dim; i++) {
			Sv[i] += v[i];

			for (unsigned int j = 0; j < dim; j++) {
				Svv[i][j] += v[i] * v[j];
			}
		}
	}

	void Fill(double v1, double v2 = 0., double v3 = 0., double v4 = 0., double v5 = 0.)
	{
		vector<double> v(5);
		v[0] = v1;
		v[1] = v2;
		v[2] = v3;
		v[3] = v4;
		v[4] = v5;

		Fill(v);
	}

	void PrintMeanAndRMS()
	{
		for (unsigned int i = 0; i < dim; i++) {
			double mu = Sv[i] / S1;
			double v = (Svv[i][i] - Sv[i]*Sv[i] / S1) / (S1 - 1.);
			double s = (v > 0.) ? sqrt(v) : 0.;
			printf("quantity %2i: mean %+.3E, std. dev. = %.3E\n", i, mu, s);
		}
	}

	void PrintCovariance()
	{
		printf("      ");
		for (unsigned int i = 0; i < dim; i++) {
			printf("   qu. %2i", i);
		}
		printf("\n");

		for (unsigned int i = 0; i < dim; i++) {

			printf("qu. %2i", i);
			for (unsigned int j = 0; j < dim; j++) {
				double c = (Svv[i][j] - Sv[i]*Sv[j] / S1) / (S1 - 1.);
				printf("   %+.3f", c);
			}
			printf("\n");
		}
	}

	void PrintCorrelation()
	{
		printf("      ");
		for (unsigned int i = 0; i < dim; i++) {
			printf("   qu. %2i", i);
		}
		printf("\n");

		for (unsigned int i = 0; i < dim; i++) {

			printf("qu. %2i", i);
			for (unsigned int j = 0; j < dim; j++) {
				double vi = (Svv[i][i] - Sv[i]*Sv[i] / S1) / (S1 - 1.);
				double vj = (Svv[j][j] - Sv[j]*Sv[j] / S1) / (S1 - 1.);
				double c = (Svv[i][j] - Sv[i]*Sv[j] / S1) / (S1 - 1.);
				double rho = c / sqrt(vi*vj);
				printf("   %+.3f", rho);
			}
			printf("\n");
		}
	}
};


