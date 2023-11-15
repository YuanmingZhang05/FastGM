#include<vector>
#include<numeric>
#include<algorithm>
#include<cmath>
class ErrorMetrics
{
public:
	ErrorMetrics() = default;
	double rmse(std::vector<double> pred, std::vector<double> target);
	double mae(std::vector<double> pred, std::vector<double> target);
};