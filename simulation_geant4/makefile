first: showers analysis dump_hits

showers: showers.cc
	g++ -O3 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		-I"${CMSSW_BASE}/src" \
		-I"/afs/cern.ch/cms/slc6_amd64_gcc481/external/fastjet/3.0.3/include" \
		-L"${CMSSW_BASE}/lib/slc6_amd64_gcc481" \
		-lTotemAnalysisTotemNtuplizer \
		showers.cc -o showers

analysis: analysis.cc
	g++ -O3 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		-I"${CMSSW_BASE}/src" \
		-I"/afs/cern.ch/cms/slc6_amd64_gcc481/external/fastjet/3.0.3/include" \
		-L"${CMSSW_BASE}/lib/slc6_amd64_gcc481" \
		-lTotemAnalysisTotemNtuplizer \
		analysis.cc -o analysis

dump_hits: dump_hits.cc
	g++ -O3 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		-I"${CMSSW_BASE}/src" \
		-I"/afs/cern.ch/cms/slc6_amd64_gcc481/external/fastjet/3.0.3/include" \
		-L"${CMSSW_BASE}/lib/slc6_amd64_gcc481" \
		-lTotemAnalysisTotemNtuplizer \
		dump_hits.cc -o dump_hits
