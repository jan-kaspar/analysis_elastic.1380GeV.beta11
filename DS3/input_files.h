#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

std::string input_ntuple_name;

void InitInputFiles()
{
	input_ntuple_name = "TotemNtuple";

	input_files.clear();

	// 2013_02_14

	/*
	std::string prefix = "rfio:///castor/cern.ch/user/j/jkaspar/reco/sr+hsx/";
	input_files.push_back(prefix + "9078_ntuple.root");
	*/

	std::string prefix = "rfio:///castor/cern.ch/totem/offline/Analysis/2013/Physics/9078/v_4.0/";
	input_files.push_back(prefix + "val9078_totem_ntuple.root");
}
