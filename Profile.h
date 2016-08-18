#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------------------------------

class Profile
{
  public:
	string name, title;

	bool uniformBinning;

	unsigned int bins;
	double from, to;
	double bin_size;

	double *binEdges;

	vector<unsigned int> S1;
	vector<double> Sv;
	vector<double> Svv;

	Profile(string _n, string _T, unsigned int _b, double _f, double _t);
	Profile(string _n, string _T, unsigned int _b, double *_be);

	signed int GetBin(double x);
	double GetBinWidth(unsigned int);

	void Fill(double x, double v);
	TH1D* GetHist(double N);
};

//----------------------------------------------------------------------------------------------------

Profile::Profile(string _n, string _T, unsigned int _b, double _f, double _t) :
		name(_n), title(_T), uniformBinning(true), bins(_b), from(_f), to(_t)
{
	for (unsigned int i = 0; i < bins; i++) {
		S1.push_back(0);
		Sv.push_back(0.);
		Svv.push_back(0.);
	}

	bin_size = (to - from) / bins;
}

//----------------------------------------------------------------------------------------------------

Profile::Profile(string _n, string _T, unsigned int _b, double *_be) :
		name(_n), title(_T), uniformBinning(false), bins(_b)
{
	for (unsigned int i = 0; i < bins; i++) {
		S1.push_back(0);
		Sv.push_back(0.);
		Svv.push_back(0.);
	}

	binEdges = new double[bins+1];
	for (unsigned int i = 0; i <= bins; i++)
		binEdges[i] = _be[i];
}

//----------------------------------------------------------------------------------------------------

signed int Profile::GetBin(double x)
{
	if (uniformBinning) {
		if (x < from || x > to)
			return -1;

		unsigned int i = (x - from) / bin_size;
		if (i >= bins)
			return -1;

		return i;
	} else {
		if (x < binEdges[0] || x > binEdges[bins])
			return -1;

		return TMath::BinarySearch(bins, binEdges, x);
	}
}

//----------------------------------------------------------------------------------------------------

void Profile::Fill(double x, double v)
{
	signed int idx = GetBin(x);

	if (idx < 0)
		return;
	
	S1[idx]++;
	Sv[idx] += v;
	Svv[idx] += v*v;
}

//----------------------------------------------------------------------------------------------------

double Profile::GetBinWidth(unsigned int idx)
{
	if (uniformBinning)
		return bin_size;
	else {
		return binEdges[idx+1] - binEdges[idx];
	}
}

//----------------------------------------------------------------------------------------------------

TH1D* Profile::GetHist(double N=1.)
{
	TH1D *h;
	if (uniformBinning)
		h = new TH1D(name.c_str(), title.c_str(), bins, from, to);
	else
		h = new TH1D(name.c_str(), title.c_str(), bins, binEdges);

	for (unsigned int i = 0; i < bins; i++) {
		double w = GetBinWidth(i);

		double M = Sv[i]/N;
		h->SetBinContent(i+1, M / w);

		double S2 = Svv[i]/N - M*M;
		double si = (S2 > 0.) ? sqrt(S2) : 0.;
		h->SetBinError(i+1, si/sqrt(N) / w);

	}

	return h;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

class Profile2D
{
  public:
	string name, title;
	
	bool uniformBinning;

	unsigned int x_bins;
	double x_from, x_to;
	double x_bin_size;
	double *x_binEdges;
	
	unsigned int y_bins;
	double y_from, y_to;
	double y_bin_size;
	double *y_binEdges;

	vector< vector<unsigned int> > S1;
	vector< vector<double> > Sv;
	vector< vector<double> > Svv;

	Profile2D(string _n, string _T, unsigned int _xb, double _xf, double _xt, unsigned int _yb, double _yf, double _yt);
	Profile2D(string _n, string _T, unsigned int _xb, double *_xbe, unsigned int _yb, double *_ybe);

	void GetBin(double x, double y, signed int &i_x, signed int &i_y);
	void Fill(double x, double y, double v);
	double GetXBinWidth(unsigned int idx);
	double GetYBinWidth(unsigned int idx);
	TH2D* GetHist(double N);
};

//----------------------------------------------------------------------------------------------------

Profile2D::Profile2D(string _n, string _T, unsigned int _xb, double _xf, double _xt, unsigned int _yb, double _yf, double _yt) :
		name(_n), title(_T), uniformBinning(true), x_bins(_xb), x_from(_xf), x_to(_xt), y_bins(_yb), y_from(_yf), y_to(_yt)
{
	for (unsigned int x_i = 0; x_i < x_bins; x_i++) {
		vector<unsigned int> s1;
		vector<double> sv;
		vector<double> svv;
		
		for (unsigned int y_i = 0; y_i < y_bins; y_i++) {
			s1.push_back(0);
			sv.push_back(0.);
			svv.push_back(0.);
		}

		S1.push_back(s1);
		Sv.push_back(sv);
		Svv.push_back(svv);
	}

	x_bin_size = (x_to - x_from) / x_bins;
	y_bin_size = (y_to - y_from) / y_bins;
}

//----------------------------------------------------------------------------------------------------

Profile2D::Profile2D(string _n, string _T, unsigned int _xb, double *_xbe, unsigned int _yb, double * /*_ybe*/) : 
	name(_n), title(_T), uniformBinning(false), x_bins(_xb), y_bins(_yb)
{
	for (unsigned int x_i = 0; x_i < x_bins; x_i++) {
		vector<unsigned int> s1;
		vector<double> sv;
		vector<double> svv;
		
		for (unsigned int y_i = 0; y_i < y_bins; y_i++) {
			s1.push_back(0);
			sv.push_back(0.);
			svv.push_back(0.);
		}

		S1.push_back(s1);
		Sv.push_back(sv);
		Svv.push_back(svv);
	}
	
	x_binEdges = new double[x_bins+1];
	for (unsigned int i = 0; i <= x_bins; i++)
		x_binEdges[i] = _xbe[i];
	
	y_binEdges = new double[y_bins+1];
	for (unsigned int i = 0; i <= y_bins; i++)
		y_binEdges[i] = _xbe[i];
}

//----------------------------------------------------------------------------------------------------

void Profile2D::GetBin(double x, double y, signed int &i_x, signed int &i_y)
{
	if (uniformBinning) {
		if (x < x_from || x > x_to)
			i_x = -1;
		else
			i_x = (x - x_from) / x_bin_size;
	
		if (y < y_from || y > y_to)
			i_y = -1;
		else
			i_y = (y - y_from) / y_bin_size;
	} else {
		if (x < x_binEdges[0] || x > x_binEdges[x_bins])
			i_x = -1;
		else
			i_x = TMath::BinarySearch(x_bins, x_binEdges, x);

		if (y < y_binEdges[0] || y > y_binEdges[y_bins])
			i_y = -1;
		else
			i_y = TMath::BinarySearch(y_bins, y_binEdges, y);
	}
}

//----------------------------------------------------------------------------------------------------

void Profile2D::Fill(double x, double y, double v)
{
	signed int x_i, y_i;
	GetBin(x, y, x_i, y_i);

	if (x_i < 0 || y_i < 0)
		return;

	S1[x_i][y_i]++;
	Sv[x_i][y_i] += v;
	Svv[x_i][y_i] += v*v;
}

//----------------------------------------------------------------------------------------------------

double Profile2D::GetXBinWidth(unsigned int idx)
{
	return (uniformBinning) ? x_bin_size : x_binEdges[idx+1] - x_binEdges[idx] ;
}

//----------------------------------------------------------------------------------------------------

double Profile2D::GetYBinWidth(unsigned int idx)
{
	return (uniformBinning) ? y_bin_size : y_binEdges[idx+1] - y_binEdges[idx] ;
}

//----------------------------------------------------------------------------------------------------

TH2D* Profile2D::GetHist(double N=1.)
{
	TH2D *h;
	if (uniformBinning)
		h = new TH2D(name.c_str(), title.c_str(), x_bins, x_from, x_to, y_bins, y_from, y_to);
	else
		h = new TH2D(name.c_str(), title.c_str(), x_bins, x_binEdges, y_bins, y_binEdges);

	for (unsigned int x_i = 0; x_i < x_bins; x_i++) {
		for (unsigned int y_i = 0; y_i < y_bins; y_i++) {
			double M = Sv[x_i][y_i]/N;
			h->SetBinContent(x_i+1, y_i+1, M);
	
			double S2 = Svv[x_i][y_i]/N - M*M;
			double si = (S2 > 0.) ? sqrt(S2) : 0.;
			h->SetBinError(x_i+1, y_i+1, si/sqrt(N));
		}
	}

	return h;
}
