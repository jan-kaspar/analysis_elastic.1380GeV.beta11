#include <map>
#include <cmath>

#include "TGraph.h"

//----------------------------------------------------------------------------------------------------

struct CumulativeHistogram
{
	std::map<double, unsigned int> data;

	void Fill(const double &v)
	{
		data[v]++;
	}

	TGraph* GetGraph(bool ascending = true)
	{
		TGraph *g = new TGraph();

		double sum = 0;

		if (ascending)
		{
			for (auto it = data.begin(); it != data.end(); ++it)
			{
				int idx = g->GetN();
				g->SetPoint(idx, it->first, sum);
				sum += it->second;
				g->SetPoint(idx + 1, it->first, sum);
			}
		} else {
			for (auto it = data.rbegin(); it != data.rend(); ++it)
			{
				int idx = g->GetN();
				g->SetPoint(idx, it->first, sum);
				sum += it->second;
				g->SetPoint(idx + 1, it->first, sum);
			}
		}

		return g;
	}
};

//----------------------------------------------------------------------------------------------------

template <class Iterator>
double KolmogorovTestLoop(Iterator beg_1, Iterator end_1, Iterator beg_2, Iterator end_2, bool ascending, TGraph *g_diff)
{
	using namespace std;

	double sum_1 = 0., sum_2 = 0.;

	Iterator it_1 = beg_1;
	Iterator it_2 = beg_2;

	double max_diff = 0.;
	while (it_1 != end_1 || it_2 != end_2)
	{
		unsigned int next = 0;

		if (it_1 == end_1)
			next = 2;

		if (it_2 == end_2)
			next = 1;

		if (next == 0)
		{
			if (it_1->first == it_2->first)
				next = 3;
			else
			{
				bool order = (it_1->first < it_2->first);
				if (!ascending)
					order = !order;
				next = (order) ? 1 : 2;
			}
		}

		double x_next = 0.; 

		if (next == 1 || next == 3)
		{
			sum_1 += it_1->second;
			x_next = it_1->first;
			++it_1;
		}

		if (next == 2 || next == 3)
		{
			sum_2 += it_2->second;
			x_next = it_2->first;
			++it_2;
		}

		double diff = sum_2 - sum_1;
		if (fabs(diff) > fabs(max_diff))
			max_diff = diff;
		
		printf("  next = %i, x_next = %.3f, diff = %.3f\n", next, x_next, diff);

		if (g_diff != NULL)
		{
			g_diff->SetPoint(g_diff->GetN(), x_next, diff);
		}
	}

	return max_diff;
}

//----------------------------------------------------------------------------------------------------

double KolmogorovTest(const CumulativeHistogram &h_1, const CumulativeHistogram &h_2, bool ascending, TGraph *g_diff = NULL)
{
	using namespace std;

	const map<double, unsigned int> &data_1 = h_1.data;
	const map<double, unsigned int> &data_2 = h_2.data;

	if (ascending)
		return KolmogorovTestLoop(data_1.begin(), data_1.end(), data_2.begin(), data_2.end(), ascending, g_diff);
	else
		return KolmogorovTestLoop(data_1.rbegin(), data_1.rend(), data_2.rbegin(), data_2.rend(), ascending, g_diff);
}
