#include "fastgm_bbm.h"
struct PoissonProcess
{
	std::vector<int> position_map;
	int id, len, n_empty_bins, n_incoming, c;
	double x;
	PoissonProcess(int i, int k);
	PoissonProcess() = default;
	void get_next_ball();
};

PoissonProcess::PoissonProcess(int i, int k):
	position_map(k), id(i), len(k), n_empty_bins(k), n_incoming(0),
	c(0), x(0.)
{
	std::iota(position_map.begin(), position_map.end(), 0);
}
void PoissonProcess::get_next_ball()
{
	float u, interval;
	int z, j;
	mmh3_generator random_engine(id, n_incoming);
	std::uniform_real_distribution<double> uniform(0., 1.);
	uint32_t hash = mmh3(id, n_incoming);
	u = uniform(random_engine);
	if (n_empty_bins < len / 10)
	{
		std::uniform_int_distribution<int> pos_dist(0, len - 1);
		j = pos_dist(random_engine);
		interval = -log(u);
		z = 1;
	}
	else
	{
		z = floor(log(u) / log(1. - (double)n_empty_bins / len)) + 1.;
		std::gamma_distribution<double> interval_dist(z, 1.);
		std::uniform_int_distribution<int> pos_dist(0, n_empty_bins - 1);
		interval = interval_dist(random_engine);
		j = pos_dist(random_engine);
	}
	c = position_map[j];
	if (j < n_empty_bins)
	{
		std::swap(position_map[j], position_map[n_empty_bins - 1]);
		n_empty_bins -= 1;
	}
	x += interval;
	n_incoming += z;
}

FastGM_BBM::FastGM_BBM():HashAlgorithm("FastGM-BBM")
{
}

Context FastGM_BBM::generate_sketch(int k, std::vector<double> v) const
{
	int n = v.size(), R_scale = 0, n_empty_bins = k, n_running = n;
	Context context(n, k);
	std::vector<bool> is_finished(n, false);
	std::vector<int> R(n);
	std::vector<PoissonProcess> generators(n);
	std::generate(generators.begin(), generators.end(), [i = 0, k]() mutable {i++; return PoissonProcess(i, k); });
	double sum_v = std::accumulate(v.begin(), v.end(), 0.);
	while (n_empty_bins > 0)
	{
		R_scale += k;
		for (int i = 0; i < n; i++)
		{
			R[i] = ceil(R_scale * v[i] / sum_v);
			while (generators[i].n_empty_bins > 0 && generators[i].n_incoming < R[i])
			{
				generators[i].get_next_ball();
				double b = generators[i].x / (k * v[i]);
				int c = generators[i].c;
				if (b < context.y[c])
				{
					context.y[c] = b;
					context.s[c] = i;
				}
				else if (context.y[c] < 0)
				{
					context.y[c] = b;
					context.s[c] = i;
					n_empty_bins -= 1;
				}
			}
		}
	}
	int argmax_y = std::max_element(context.y.begin(), context.y.end()) - context.y.begin();
	while (n_running)
	{
		R_scale += k;
		for (int i = 0; i < n; i++)
		{
			if (is_finished[i])
				continue;
			R[i] = ceil(R_scale * v[i] / sum_v);
			while (generators[i].n_empty_bins > 0 && generators[i].n_incoming < R[i])
			{
				generators[i].get_next_ball();
				double b = generators[i].x / (k * v[i]);
				int c = generators[i].c;
				if (b > context.y[argmax_y])
				{
					is_finished[i] = true;
					n_running -= 1;
					break;
				}
				if (b < context.y[c])
				{
					context.y[c] = b;
					context.s[c] = i;
					if (c == argmax_y)
					{
						argmax_y = std::max_element(context.y.begin(), context.y.end()) - context.y.begin();
					}
				}
			}
		}
	}
	return context;
}
