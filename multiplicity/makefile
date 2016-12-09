all: analyze

analyze: analyze.cc
	g++ -O3 -g -Wall -Wextra --std=c++11\
			-I../ \
			`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
			-I$(CMSSW_BASE)/src \
			-I/afs/cern.ch/cms/slc6_amd64_gcc481/external/fastjet/3.0.3/include \
			-L$(CMSSW_BASE)/lib/slc6_amd64_gcc481 \
			-lTotemAnalysisTotemNtuplizer \
			analyze.cc -o analyze
