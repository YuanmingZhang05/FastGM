#include<cstdint>
uint32_t mmh3(uint32_t key, int seed);
double normalize(uint32_t hash);
class mmh3_generator
{
public:
	uint32_t key_, seed_;
	using result_type = uint32_t;
	mmh3_generator(uint32_t key, uint32_t seed);
	static uint32_t min() { return 0; }
	static uint32_t max() { return (uint32_t)-1; }
	uint32_t operator()();
};