#pragma once

#include <vector>

class DesignVariable
{
public:
	DesignVariable(int length, double minVal, double maxVal) : m_length{ length }, m_minVal{ minVal }, m_maxVal{ maxVal } {};
	DesignVariable(double precision, double minVal, double maxVal) : m_precision{ precision }, m_minVal{ minVal }, m_maxVal{ maxVal } {};
	
	int GetLength() { return m_length; }
	double GetPrecision() { return m_precision; }
	double GetMinValue() { return m_minVal; }
	double GetMaxValue() { return m_maxVal; }

private:

	int m_length = -INT_MAX;
	double m_precision = -DBL_MAX;
	double m_minVal = -DBL_MAX;
	double m_maxVal = DBL_MAX;
};

