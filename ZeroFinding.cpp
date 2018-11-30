/* ZeroFinding.cpp 
Description:
	* 




*/

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>
#include "ZeroFinding.hpp"


#pragma region Constructors/Destructor:
ZeroFinding::ZeroFinding() : Results()
{

}
ZeroFinding::~ZeroFinding()
{

}
#pragma endregion
#pragma region Class Methods
double ZeroFinding::FindZero_Bisect(FType func, double a, double b, double tol_int, double tol_approx, std::vector<std::tuple<unsigned, double, double>> *tracker)
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Description:
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Return approximate zero to one-dimensional scalar function func using bisection method.
	if (tracker == nullptr)
	{
		this->Results.clear();
	}
	double x_l, x_r, x_m;
	unsigned n = 1;
	x_l = a;
	x_r = b;
	// Use algorithm using stopping criterion:
	while (std::max(std::abs(func(x_l)), std::abs(func(x_r))) > tol_approx || (x_r - x_l) > tol_int)
	{
		x_m = (x_l + x_r) / 2.0;
		// Set active interval based upon sign of function:
		if (func(x_l) * func(x_m) < 0)
		{
			x_r = x_m;
		}
		else
		{
			x_l = x_m;
		}
		if (tracker != nullptr)
		{
			// Append the step number, x_m and function approximation to the tracking vector if passed:
			tracker->push_back(std::make_tuple(x_m, func(x_m), n++));
		}
		else
		{
			this->Results.push_back(std::make_tuple(n++, x_m, func(x_m)));
		}
	}
	// Return the approximate zero within given tolerances:
	return x_m;
}
double ZeroFinding::FindZero_Newtons(FType func, FType firstDeriv, double x_0, double tol_consec, double tol_approx, std::vector<std::tuple<unsigned, double, double>> *tracker)
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Description:
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Return approximate zero to one-dimensional scalar function func using Newton's method
	// (uses true first derivative of underlying function)
	if (tracker == nullptr)
	{
		this->Results.clear();
	}
	double x_new = x_0;
	double x_old = x_0 - 1;
	unsigned n = 1;
	while (std::abs(func(x_new)) > tol_approx || (std::abs(x_new - x_old) > tol_consec))
	{
		// Throw exception if the first derivative is zero at evaluation, since cannot divide by zero:
		if (!(firstDeriv(x_new)))
		{
			throw std::exception("Cannot divide by zero.");
		}
		x_old = x_new;
		x_new = x_old - (func(x_old) / firstDeriv(x_old));
		if (tracker != nullptr)
		{
			// Append the step number, x_m and function approximation to the tracking vector if passed:
			tracker->push_back(std::make_tuple(x_new, func(x_new), n++));
		}
		else
		{
			this->Results.push_back(std::make_tuple(n++, x_new, func(x_new)));
		}
	}
	// Return the approximate zero:
	return x_new;
}

double ZeroFinding::FindZero_Secant(FType func, double x_0, double x_minus_1, double tol_consec, double tol_approx, std::vector<std::tuple<unsigned, double, double>> *tracker) 
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Description:
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Return approximate zero to one-dimensional scalar function func using the secant method 
	// (uses approximation of first order derivative using finite difference method)
	if (tracker == nullptr)
	{
		this->Results.clear();
	}
	double x_new = x_0;
	double x_old = x_minus_1;
	double x_oldest;
	unsigned n = 1;
	while (std::abs(func(x_new)) > tol_approx || std::abs(x_new - x_old) > tol_consec)
	{
		x_oldest = x_old;
		x_old = x_new;
		// Throw exception if func(x_old) - func(x_oldest) is zero:
		if (!(func(x_old) - func(x_oldest)))
		{
			throw std::exception("Cannot divide by zero.");
		}
		x_new = x_old - func(x_old) * (x_old - x_oldest) / (func(x_old) - func(x_oldest));
		if (tracker != nullptr)
		{
			// Append the step number, x_m and function approximation to the tracking vector if passed:
			tracker->push_back(std::make_tuple(x_new, func(x_new), n++));
		}
		else
		{
			this->Results.push_back(std::make_tuple(n++, x_new, func(x_new)));
		}
	}
	return x_new;
}
void ZeroFinding::PrintResults(unsigned precision, const std::string &FunctionName)
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Description:
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Print the stored results of the most recent zero finding method.
	if (this->Results.size() > 0)
	{
		if (FunctionName != "")
		{
			// Print function name:
			std::cout << "F(X) = " << FunctionName << std::endl;
		}
		std::cout << "N:" << std::setw(precision) << "x_n:" << std::setw(precision) << "F(x):" << std::endl;
		for (auto iter = this->Results.begin(); iter != this->Results.end(); iter++)
		{
			std::cout << std::setprecision(0) << std::get<0>(*iter) << "\t" << std::setprecision(precision) << std::fixed << std::get<1>(*iter) << "\t" << std::get<2>(*iter) << std::endl;
		}
	}
}
#pragma endregion
#pragma region Overloaded Operators
ZeroFinding& ZeroFinding::operator=(const ZeroFinding &in)
{
	if (this != &in)
	{
		this->Results = in.Results;
	}
	return *this;
}
#pragma endregion


