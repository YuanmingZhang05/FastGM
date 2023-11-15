#include"estimators.h"
#include<iostream>
double Estimator::estimate_lambda(std::vector<double> s) const
{
	double G_m = std::accumulate(s.begin(), s.end(), 0.);
	return (s.size() - 1.) / G_m;
}
double Estimator::lambda(std::vector<double> v) const
{
	return std::accumulate(v.begin(), v.end(), 0.);
}
double Estimator::estimate_J_p(std::vector<int> s1, std::vector<int> s2) const
{
	return std::inner_product(s1.begin(), s1.end(), s2.begin(), 0.,
		std::plus<>(), [](int a, int b) -> double {return a == b ? 1. : 0.; }) / s1.size();
}
double Estimator::J_p(std::vector<double> v1, std::vector<double> v2) const
{
	double j_p = 0;
	for (int i = 0; i < v1.size(); i++)
	{
		double t = 0.;
		for (int j = 0; j < v1.size(); j++)
		{
			t += std::max(v1[j] / v1[i], v2[j] / v2[i]);
		}
		j_p += 1. / t;
	}
	return j_p;

}