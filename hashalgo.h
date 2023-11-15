#pragma once
#include<vector>
#include<string>
#include"mmh3.h"
class Context
{
public:
	std::vector<double> y;
	std::vector<int> s;
	Context(int n, int k);
	Context() = default;
};
class HashAlgorithm
{
public:
	HashAlgorithm(std::string name);
	std::string algorithm_name;
	virtual Context generate_sketch(int k, std::vector<double> v) const = 0;
};