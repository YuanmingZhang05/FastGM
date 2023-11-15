#include<vector>
#include<numeric>
#include<random>
#include<cmath>
#include<algorithm>
#include<cstdlib>
#include"hashalgo.h"
class FastGM_BBM : public HashAlgorithm
{
public:
	FastGM_BBM();
	Context generate_sketch(int k, std::vector<double> v) const;
};