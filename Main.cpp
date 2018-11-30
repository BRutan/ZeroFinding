/* Main.cpp
Description:
	* 




*/
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <iostream>
#include "ZeroFinding.hpp"

int main()
{
	////////////////////////////////////
	// Problem 5.1:
	////////////////////////////////////
	double S_0, K, r, q, t;
	double SQRT_PI_1 = 1 / std::sqrt(2.0 * M_PI);
	double fix = 1 / (2 * std::sqrt(2.0));
	double Price = 2.5;
	S_0 = 30;
	K = 30;
	r = .03;
	q = .01;
	t = .5;
	FType d_1 = [&S_0, &K, &r, &q, &t](double iV) { return (std::log(S_0 / K) + ((r - q + (iV * iV / 2.0)) * t)) / (iV * std::sqrt(t)); };
	FType d_2 = [&d_1, &S_0, &K, &r, &q, &t](double iV) { return d_1(iV) - iV * std::sqrt(t); };
	FType Norm_CDF = [&fix](double Z) { return 0.5 * erfc(-Z * M_SQRT1_2); };
	FType Norm_PDF = [&SQRT_PI_1](double Z) { return SQRT_PI_1 * std::exp(-Z * Z * .5); };
	FType CallOptionPrice = [&d_1, &d_2, &Norm_CDF, &S_0, &K, &r, &q, &t, &Price](double iV) { return (S_0 * std::exp(-q * t) * Norm_CDF(d_1(iV))) - (K * std::exp(-r * t) * Norm_CDF(d_2(iV))) - Price; };
	FType Vega = [&d_1, &d_2, &Norm_PDF, &S_0, &K, &r, &q, &t](double iV) { return S_0 * std::exp(-q * t) * std::sqrt(t) * Norm_PDF(d_1(iV)); };
	/////////////////
	// Compute the implied volatility for this call option using Newton's Method:
	/////////////////
	double tol_approx = 10E-6;
	double tol_consec = 10E-6;
	ZeroFinding finder;
	finder.FindZero_Newtons(CallOptionPrice, Vega, 0.5, tol_consec, tol_approx);
	finder.PrintResults(16);
	////////////////////////////////////
	// Problem 5.3:
	////////////////////////////////////
	// Compute using secant method, newton's method and bisection method for different option:
	double result;
	S_0 = 40.0;
	K = 40.0;
	q = .01;
	r = .025;
	t = 5.0f / 12.0f;
	Price = 2.75;
	// Compute using bisection method on interval [0.0001,1]:
	result = finder.FindZero_Bisect(CallOptionPrice,0.0001, 1, tol_approx, tol_consec);
	std::cout << "Bisection method: " << std::endl;
	finder.PrintResults(6);
	// Compute using secant method:
	result = finder.FindZero_Secant(CallOptionPrice, 0.5, 0.5 - 1, tol_consec, tol_approx);
	std::cout << "Secant method: " << std::endl;
	finder.PrintResults(6);
	// Compute using Newton's method:
	result = finder.FindZero_Newtons(CallOptionPrice, Vega, 0.5, tol_consec, tol_approx);
	std::cout << "Newton's method: " << std::endl;
	finder.PrintResults(6);
	////////////////////////////////////
	// Problem 5.4
	////////////////////////////////////
	// Find yield of three year semiannual coupon bond with 4% coupon rate, price 101:
	double step = 0.5;
	double c = .04;
	t = 3.0;
	Price = 1.01;
	FType PresentValue = [&t, &step, &Price, &c](double ytm) 
	{ 
		double output = 0.0;
		double start = step;
		// Generate the sum of the discounted cash flow stream:
		while (start < t)
		{
			output += c * step * std::exp(-ytm * start);
			start += step;
		}
		output += (c * step + 1) * std::exp(-ytm * start);
		return output - Price; 
	};
	FType duration = [&t, &step, &Price, &c](double ytm)
	{
		double output = 0.0;
		double start = step;
		// Generate the sum of the discounted cash flow_i * maturity_i stream:
		while (start < t)
		{
			output += -c * step * start * std::exp(-ytm * start);
			start += step;
		}
		output += -(c * step + 1) * start * std::exp(-ytm * start);
		return output;
	};
	FType convexity = [&t, &step, &Price, &c](double ytm)
	{
		double output = 0.0;
		double start = step;
		// Generate the sum of the discounted cash flow_i * maturity_i^2 stream:
		while (start < t)
		{
			output += c * step * start * start * std::exp(-ytm * start);
			start += step;
		}
		output += (c * step + 1) * start * start * std::exp(-ytm * start);
		return output;
	};
	// Compute the ytm, modified duration and convexity of the bond using Newton's method:
	double ytm_output = finder.FindZero_Newtons(PresentValue, duration, .04, tol_consec, tol_approx);
	finder.PrintResults(6, "Yield to Maturity");
	std::cout << "Modified Duration: " << -duration(ytm_output) / Price << std::endl;
	std::cout << "Convexity: " << convexity(ytm_output) << std::endl;

	system("pause");

	return 0;
}