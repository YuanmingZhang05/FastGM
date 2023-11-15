#include"metrics.h"
double ErrorMetrics::rmse(std::vector<double> pred, std::vector<double> target)
{
	return sqrt(std::inner_product(pred.begin(), pred.end(), target.begin(), 0.,
		std::plus<>(), [](double pred, double target) -> double {return pow(pred - target, 2); }) / pred.size());
}
double ErrorMetrics::mae(std::vector<double> pred, std::vector<double> target)
{
	return std::inner_product(pred.begin(), pred.end(), target.begin(), 0.,
		std::plus<>(), [](double pred, double target) -> double {return abs(pred - target); }) / pred.size();
}
