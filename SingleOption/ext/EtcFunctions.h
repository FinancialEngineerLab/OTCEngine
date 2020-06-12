#ifndef ETCFUNCTIONS_H
#define ETCFUNCTIONS_H
#include <vector>

class EtcFunctions{
public:
	EtcFunctions(){};
	double intp1d(std::vector<double> v, std::vector<double> axis, double target) const;
	
};
#endif