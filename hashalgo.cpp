#include"hashalgo.h"
Context::Context(int n, int k) : y(k, -1.), s(k)
{
}

HashAlgorithm::HashAlgorithm(std::string name) : algorithm_name(name)
{
}
