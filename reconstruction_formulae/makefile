first: test_formulae test_formulae_graph

test_formulae: test_formulae.cc formulae.h ../common_definitions.h ../common_parameters.h ../common_algorithms.h
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall test_formulae.cc -o test_formulae

test_formulae_graph: test_formulae_graph.cc formulae.h ../common_definitions.h ../common_parameters.h ../common_algorithms.h
	g++ `root-config --libs` -lMinuit `root-config --cflags` -O3 -Wall test_formulae_graph.cc -o test_formulae_graph
