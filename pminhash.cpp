#include"pminhash.h"
#include<iostream>
PMinHash::PMinHash() : HashAlgorithm("PMinHash")
{
}
Context PMinHash::generate_sketch(int k, std::vector<double> v) const
{
	int n = v.size();
	Context context(n, k);
	std::uniform_real_distribution<double> uniform;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < k; j++)
		{
			mmh3_generator random_engine(i, j);
			double u = uniform(random_engine), x;
			x = -log(u) / v[i];
			if (x < context.y[j] || context.y[j] < 0)
			{
				context.y[j] = x;
				context.s[j] = i;
			}

		}
	}
	return context;
}