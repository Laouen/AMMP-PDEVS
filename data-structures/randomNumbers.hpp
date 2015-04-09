#ifndef BOOST_SIMULATION_PDEVS_RANDOMNUMBERS_H
#define BOOST_SIMULATION_PDEVS_RANDOMNUMBERS_H

#include <random>
using namespace std;

template<class NumbType>
class IntegerRandom {

public:

	IntegerRandom() = default;
	IntegerRandom(std::mt19937::result_type s) : generator(s) {};
	
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
class RealRandom {

public:

	RealRandom() = default;
	RealRandom(std::mt19937::result_type s) : generator(s) {};
	
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