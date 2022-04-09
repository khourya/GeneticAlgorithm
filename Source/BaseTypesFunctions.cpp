#include <iostream>

#include "BaseTypesFunctions.h"

bool BaseFunctions::WeightedCoinFlip(double p)
{
	if (p == 0.)
		return false;

	if (p == 1.)
		return true;

	double coinFlipProbability = static_cast<double>(rand()) / static_cast<double> (RAND_MAX);
	return coinFlipProbability <= p ? true : false;
}
