#include<vector>
#include<numeric>
#include<random>
#include<cmath>
#include<algorithm>
#include<cstdlib>
#include"hashalgo.h"
class FastGM : public HashAlgorithm
{
public:
	FastGM();
	Context generate_sketch(int k, std::vector<double> v) const;
};