#include<vector>
#include<cmath>
#include<algorithm>
#include<random>
#include"hashalgo.h"
#include<cstdlib>
class PMinHash : public HashAlgorithm
{
public:
	PMinHash();
	Context generate_sketch(int k, std::vector<double> v) const;
};