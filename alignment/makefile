first: test1 test2 fit_th_y_diff match_th_y_scale_fcn match_th_y_crop match_th_y_scale_hist kolmogorov_validation \
		simulation compile_results

test1: test1.cc
	g++ -O3 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		test1.cc -o test1

test2: test2.cc
	g++ -O3 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		test2.cc -o test2

fit_th_y_diff: fit_th_y_diff.cc
	g++ -O3 --std=c++11 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		fit_th_y_diff.cc -o fit_th_y_diff

match_th_y_scale_fcn: match_th_y_scale_fcn.cc
	g++ -O3 --std=c++11 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		match_th_y_scale_fcn.cc -o match_th_y_scale_fcn

match_th_y_crop: match_th_y_crop.cc kolmogorov.h
	g++ -O3 --std=c++11 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		match_th_y_crop.cc -o match_th_y_crop

match_th_y_scale_hist: match_th_y_scale_hist.cc
	g++ -O3 --std=c++11 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		match_th_y_scale_hist.cc -o match_th_y_scale_hist

kolmogorov_validation: kolmogorov_validation.cc kolmogorov.h
	g++ -O3 --std=c++11 -Wall -Wextra\
		`root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` \
		kolmogorov_validation.cc -o kolmogorov_validation

simulation: simulation.cc
	g++ -O3 --std=c++11 -Wall -Wextra\
		`root-config --libs` `root-config --cflags` \
		simulation.cc -o simulation

compile_results: compile_results.cc
	g++ -O3 --std=c++11 -Wall -Wextra\
		`root-config --libs` `root-config --cflags` \
		compile_results.cc -o compile_results
