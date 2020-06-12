#include <stdexcept>
#include <vector> 
#include "EtcFunctions.h"



double EtcFunctions::intp1d(std::vector<double> v, std::vector<double> axis, double t) const
{
	//check size
	if(v.size()!=axis.size()) throw std::logic_error("size mismatched");

	double result;

	if (t<=axis.front())
	{
		t=axis.front();
	}
    else if (t>=axis.back())
	{
		t=axis.back();
	}
	
	std::vector<double>::iterator iter_v=v.begin()+1;
	std::vector<double>::iterator iter_x=axis.begin()+1;

	for(;iter_x != axis.end();iter_x++)
	{
	   if (t <= *iter_x)
		{
            result = (*(iter_v) * (t - *(iter_x-1)) + *(iter_v-1) * (*iter_x - t)) / (*iter_x- *(iter_x-1));
            return result;
		}
	   iter_v++;
	}

	
	throw std::logic_error("interpolation fails");
	return 1;

}

