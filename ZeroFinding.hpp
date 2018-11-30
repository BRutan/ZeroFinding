/* ZeroFinding.hpp
Description:
	* Class performs three zero finding algorithms for given functions, using Bisection, Newton's and Secant method.
Class Members:
	* Vector<tuple<unsigned, double, double>> Results: Store results of most recent operation. Will not be filled if a similar tracker vector is passed to function.
Class Methods:
	* double FindZero_Bisect(): Find zero using the bisection method.
	* 





*/


#ifndef ZEROFINDING_HPP
#define ZEROFINDING_HPP

#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using FType = std::function<double(double)>;

class ZeroFinding
{
private:
	std::vector<std::tuple<unsigned, double, double>> Results;
public:
	/////////////////////////////
	// Constructors/Destructor:
	/////////////////////////////
	ZeroFinding();
	virtual ~ZeroFinding();
	/////////////////////////////
	// Class Methods:
	/////////////////////////////
	double FindZero_Bisect(FType func, double a, double b, double tol_int, double tol_approx, std::vector<std::tuple<unsigned, double, double>> *tracker = nullptr);
	double FindZero_Newtons(FType func, FType firstDeriv, double x_0, double tol_consec, double tol_approx, std::vector<std::tuple<unsigned, double, double>> *tracker = nullptr);
	double FindZero_Secant(FType func, double x_0, double x_minus_1, double tol_consec, double tol_approx, std::vector<std::tuple<unsigned, double, double>> *tracker = nullptr);
	void PrintResults(unsigned precision, const std::string &FunctionName = "");
	/////////////////////////////
	// Overloaded Operators:
	/////////////////////////////
	ZeroFinding& operator=(const ZeroFinding&);
};


#endif