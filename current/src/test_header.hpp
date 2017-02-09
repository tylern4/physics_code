// Auto Generated header file
// Made from makeHeader.cpp

#ifndef AUTO_GEN_HEADER_H
#define AUTO_GEN_HEADER_H

double gaussian_test(double x) {
	static const float inv_sqrt_2pi = 0.3989422804014327;
	a = 1.000000;
	m = 20.000000;
	s = 3.000000;
	double p = (x - m) / s;
	return inv_sqrt_2pi / s * a * std::exp(-0.5f * p * p);
}

#endif

