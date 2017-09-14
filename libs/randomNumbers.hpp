#ifndef BOOST_SIMULATION_PDEVS_RANDOMNUMBERS_H
#define BOOST_SIMULATION_PDEVS_RANDOMNUMBERS_H

#include <random>
using namespace std;

template<class NumbType>
class IntegerRandom_t {

public:

	IntegerRandom_t() = default;
	IntegerRandom_t(std::mt19937::result_type s) : generator(s) {};
	
	NumbType drawNumber(NumbType a, NumbType b) {
		return uniform_int_distribution<NumbType>(a, b)(generator);
	}

	void seed(std::mt19937::result_type s) {
		generator.seed(s);
	}

private:
	mt19937 generator;
};

template<class NumbType>
class RealRandom_t {

public:

	RealRandom_t() = default;
	RealRandom_t(std::mt19937::result_type s) : generator(s) {};
	
	NumbType drawNumber(NumbType a, NumbType b) {
		return uniform_real_distribution<NumbType>(a, b)(generator);
	}

	void seed(std::mt19937::result_type s) {
		generator.seed(s);
	}

private:
	mt19937 generator;
};

#endif // BOOST_SIMULATION_PDEVS_RANDOMNUMBERS_H