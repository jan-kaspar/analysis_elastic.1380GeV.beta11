#include <cstdio>
#include <cstring>
#include <string>
#include <map>

#include "../Stat.h"


using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: compile_results <option> <option> ...\n");
	printf("OPTIONS:\n");
	printf("    --inputDir\n");
	printf("    --inputFile\n");
	printf("    --from\n");
	printf("    --to\n");
}

//----------------------------------------------------------------------------------------------------

void ProcessOneFile(const string &fn, map<string, Stat> &stat)
{
	FILE *f = fopen(fn.c_str(), "r");
	if (f == NULL)
	{
		printf("ERROR: Can't open file '%s', skipping.\n", fn.c_str());
		return;
	}

	while (!feof(f))
	{
		char buf[200];
		char *ret = fgets(buf, 200, f);

		if (ret == NULL)
			continue;

		char *p = strstr(buf, ",");
		*p = 0;

		double val = atof(p+1);

		//printf("%s -> %E\n", buf, val);

		auto it = stat.find(buf);
		if (it == stat.end())
			it = stat.insert({buf, Stat(1)}).first;
		
		it->second.Fill(val);
	}

	fclose(f);
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string inputDir = "./";
	string inputFile = "bla";
	unsigned int subDirFrom = 1;
	unsigned int subDirTo = 1;

	// parse command line
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
		{
			PrintUsage();
			return 0;
		}

		if (strcmp(argv[i], "--inputDir") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--inputDir' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			inputDir = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "--inputFile") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--inputFile' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			inputFile = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "--from") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--from' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			subDirFrom = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "--to") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("ERROR: option '--to' expects a parameter.\n");
				PrintUsage();
				return 1;
			}

			subDirTo = atoi(argv[++i]);
			continue;
		}

		printf("ERROR: unknown option '%s'\n", argv[i]);
		PrintUsage();
		return 1;
	}

	// prepare data containers
	map<string, Stat> stat;

	// process data
	for (unsigned int sdi = subDirFrom; sdi <= subDirTo; sdi++)
	{
		char buf[300];
		sprintf(buf, "%s/%i/%s", inputDir.c_str(), sdi, inputFile.c_str());

		ProcessOneFile(buf, stat);
	}

	// print results
	for (const auto &p : stat)
	{
		printf("%s : N = %.0f, mean = %.3E +- %.3E, std. dev. = %.3E +- %.3E\n", p.first.c_str(),
			p.second.GetEntries(), p.second.GetMean(0), p.second.GetMeanUnc(0), p.second.GetStdDev(0), p.second.GetStdDevUnc(0));
	}

	return 0;
}
