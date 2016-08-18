#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

std::string input_ntuple_name;

void InitInputFiles()
{
	input_ntuple_name = "TotemNtuple";

	input_files.clear();

	// 2013_02_11

	std::string prefix = "rfio:///castor/cern.ch/user/j/jkaspar/reco/sr+hsx/";
	input_files.push_back(prefix + "9008_ntuple.root");
	input_files.push_back(prefix + "9009_ntuple.root");
	input_files.push_back(prefix + "9010_ntuple.root");
	
	//std::string prefix = "rfio:///castor/cern.ch/totem/offline/Analysis/2013/Physics/9009_9010/v_4.0/";
}
