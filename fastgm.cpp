#include"fastgm.h"
#include<iostream>
FastGM::FastGM() : HashAlgorithm("FastGM")
{
}
Context FastGM::generate_sketch(int k, std::vector<double> v) const
{
	int n = v.size(), n_empty_bins = k, R_scale = 0;
	double sum_v = std::accumulate(v.begin(), v.end(), 0.);
	Context context(n, k);
	std::vector<int> pos_map_proto(k), R(n), n_generated(n, 0);
	std::iota(pos_map_proto.begin(), pos_map_proto.end(), 0);
	std::vector<std::vector<int>> pos_map(n, pos_map_proto);
	std::vector<double> b(n, 0.);
	std::uniform_real_distribution<double> uniform;
	while (n_empty_bins > 0)
	{
		R_scale += k;
		for (int i = 0; i < n; i++)
		{
			R[i] = ceil(R_scale * v[i] / sum_v);
			while (n_generated[i] < R[i])
			{
				n_generated[i] += 1;
				mmh3_generator random_engine(i, n_generated[i]);
				std::uniform_int_distribution<int> pos_dist(n_generated[i] - 1, k - 1);
				double u = uniform(random_engine);
				int j = pos_dist(random_engine), pos;
				b[i] -= 1. / v[i] * log(u) / (k - n_generated[i] + 1);
				pos = pos_map[i][j];
				std::swap(pos_map[i][j], pos_map[i][n_generated[i] - 1]);
				if (context.y[pos] < 0.)
				{
					context.y[pos] = b[i];
					context.s[pos] = i;
					if((n_empty_bins -= 1) == 0)
						break;
				}
				else if (context.y[pos] > b[i])
				{
					context.y[pos] = b[i];
					context.s[pos] = i;
				}
			}
			if (n_empty_bins == 0)
				break;
		}
	}
	int argmax_y = std::max_element(context.y.begin(), context.y.end()) - context.y.begin(), n_running = n;
	std::vector<bool> running(n, true);
	while (n_running)
	{
		R_scale += k;
		for (int i = 0; i < n; i++)
		{
			if (!running[i])continue;
			R[i] = ceil(R_scale * v[i] / sum_v);
			while (n_generated[i] < R[i])
			{
				n_generated[i] += 1;
				mmh3_generator random_engine(i, n_generated[i]);
				std::uniform_int_distribution<int> pos_dist(n_generated[i] - 1, k - 1);
				double u = uniform(random_engine);
				int j = pos_dist(random_engine), pos;
				b[i] -= 1. / v[i] * log(u) / (k - n_generated[i] + 1);
				pos = pos_map[i][j];
				std::swap(pos_map[i][j], pos_map[i][n_generated[i] - 1]);
				if (b[i] > context.y[argmax_y])
				{
					running[i] = false;
					n_running -= 1;
					break;
				}
				else if (context.y[pos] > b[i])
				{
					context.y[pos] = b[i];
					context.s[pos] = i;
					if (pos == argmax_y)
					{
						argmax_y = std::max_element(context.y.begin(), context.y.end()) - context.y.begin();
					}
				}
			}
		}
	}
	return context;
}