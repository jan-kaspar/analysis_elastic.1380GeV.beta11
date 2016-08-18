#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

void InitInputFiles()
{
	input_files.clear();

	// 2013_02_12
	std::string prefix = "rfio:///castor/cern.ch/user/j/jkaspar/reco/sr+hsx/";

	input_files.push_back(prefix + "9047_ntuple.root");
	input_files.push_back(prefix + "9048_ntuple.root");
	input_files.push_back(prefix + "9049_ntuple.root");
	input_files.push_back(prefix + "9050_ntuple.root");
}
