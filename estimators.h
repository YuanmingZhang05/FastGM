#include<vector>
#include<algorithm>
#include<numeric>
class Estimator
{
public:
	Estimator() = default;
	double estimate_lambda(std::vector<double> s) const;
	double estimate_J_p(std::vector<int> s1, std::vector<int> s2) const;
	double lambda(std::vector<double> v) const;
	double J_p(std::vector<double> v1, std::vector<double> v2) const;
};