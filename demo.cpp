#include<iostream>
#include<string>
#include<functional>
#include"fastgm.h"
#include"pminhash.h"
#include"fastgm_bbm.h"
#include"estimators.h"
#include"metrics.h"
std::vector<double> random_vector(int n)
{
	std::vector<double> v(n);
	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_real_distribution<double> uniform_real(0., 1.);
	std::generate(v.begin(), v.end(), std::bind(uniform_real, random_engine));
	return v;
}
struct Metrics
{
public:
	double time_consumed_sec;
	double rmse;
	double mae;
	std::string algo;
	void print()
	{
		std::cout << algo << std::endl;
		std::cout << "Avg. running time: " << time_consumed_sec << std::endl;
		std::cout << "RMSE of Jp estimation: " << rmse << std::endl;
		std::cout << "MAE of ¦Ë estimation: " << mae << std::endl;
	}
};
void Evaluate(HashAlgorithm *algo, int k, std::vector<std::vector<double>> data)
{
	Metrics metrics;
	int b = data.size();
	std::vector<Context> sketches; 
	std::vector<double> Jp_pred, Jp, lambda_pred, lambda;
	Estimator estimator;
	ErrorMetrics error_metrics;
	std::transform(data.begin(), data.end(), std::inserter(lambda, lambda.begin()), [=](std::vector<double> v) -> double {return estimator.lambda(v); });
	time_t t_begin = clock();
	std::transform(data.begin(), data.end(), std::inserter(sketches, sketches.begin()),
		[=](std::vector<double> v) -> Context {return algo->generate_sketch(k, v); });
	time_t t_end = clock();
	std::transform(sketches.begin(), sketches.end(), std::inserter(lambda_pred, lambda_pred.begin()), [=](Context ctx) -> double {return estimator.estimate_lambda(ctx.y); });
	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<int> random_id(0, b - 1);
	int n_samples = 10000; 
	for (int sample = 0; sample < n_samples; sample++)
	{
		int i = random_id(random_engine), j = random_id(random_engine);
		Jp_pred.push_back(estimator.estimate_J_p(sketches[i].s, sketches[j].s));
		Jp.push_back(estimator.J_p(data[i], data[j]));
	}
	/*
	for (int i = 0; i < b; i++)
	{
		for (int j = 0; j < i; j++)
		{
			Jp_pred.push_back(estimator.estimate_J_p(sketches[i].s, sketches[j].s));
			Jp.push_back(estimator.J_p(data[i], data[j]));
		}
	}
	*/
	metrics.algo = algo->algorithm_name;
	metrics.rmse = error_metrics.rmse(Jp_pred, Jp);
	metrics.mae = error_metrics.mae(lambda_pred, lambda);
	metrics.time_consumed_sec = (double)(t_end - t_begin) / CLOCKS_PER_SEC / b;
	metrics.print();
}
int main(void)
{
	int n = 1000, k = 1 << 12, b = 1000;
	FastGM fastgm;
	PMinHash pminhash;
	FastGM_BBM fastgm_bbm;
	std::vector<std::vector<double>> v_matrix(b, std::vector<double>(n));
	std::generate(v_matrix.begin(), v_matrix.end(), std::bind(random_vector, n));
	Evaluate(&fastgm, k, v_matrix);
	Evaluate(&fastgm_bbm, k, v_matrix);
	Evaluate(&pminhash, k, v_matrix);
}