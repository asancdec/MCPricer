// (C++) Monte Carlo Option Pricer with Euler - Maruyama Discretization
// EuropeanOption.cpp
// Álvaro Sánchez de Carlos
// Description: this file contains the source code for the EuropeanOption class

#include <iostream>
#include <string>
#include <vector>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/constants/constants.hpp>
#include "EuropeanOption.hpp"

// Function to calculate d1 in the Black-Scholes formula
double EuropeanOption::D1() const
{
	return (log(m_S / m_K) + (m_b + m_sigma * m_sigma / 2) * m_T) / (m_sigma * std::sqrt(m_T));
}

// Function to calculate d2 in the Black-Scholes formula
double EuropeanOption::D2(const double& d1) const
{
	return d1 - m_sigma * std::sqrt(m_T);
}

// Cumulative distribution function for the standard normal distribution
double EuropeanOption::N(double x) const
{
	return boost::math::cdf(boost::math::normal(), x);
}

// Probability density function for the standard normal distribution
double EuropeanOption::N_prime(double x) const
{
	return 1 / sqrt(2 * boost::math::constants::pi<double>()) * exp(-(x * x / 2));
}

// Constructor for EuropeanOption class
EuropeanOption::EuropeanOption(const std::string& type, const double& T, const double& K, const double& S, const double& r,
	const double& sigma, const int& id, const double& b) :
	m_type(type),
	m_T(T),
	m_K(K),
	m_S(S),
	m_r(r),
	m_sigma(sigma),
	m_id(id),
	m_b(std::isnan(b) ? r : b)

{}

// Copy constructor for EuropeanOption class
EuropeanOption::EuropeanOption(const EuropeanOption& source) :
	m_type(source.m_type),
	m_T(source.m_T),
	m_K(source.m_K),
	m_S(source.m_S),
	m_r(source.m_r),
	m_sigma(source.m_sigma),
	m_id(source.m_id),
	m_b(source.m_b)

{}

// Assignment operator for EuropeanOption class
EuropeanOption& EuropeanOption::operator=(const EuropeanOption& source)
{
	// Check for self assignment
	// Compare the address of the object with 
	// the address of its source object 
	if (this == &source)
		return *this;

	m_type = source.m_type;
	m_T = source.m_T;
	m_K = source.m_K;
	m_S = source.m_S;
	m_r = source.m_r;
	m_sigma = source.m_sigma;
	m_id = source.m_id;
	m_b = source.m_b;

	// Return the current object (the object pointed by this)
	return *this;
}

// Convert option parameters to a vector of strings
std::vector<std::string> EuropeanOption::ConvertToVectorString() const
{
	return
	{
		std::to_string(m_id),
		m_type,
		std::to_string(m_T),
		std::to_string(m_K),
		std::to_string(m_S),
		std::to_string(m_r),
		std::to_string(m_sigma),
		std::to_string(m_b)
	};
}

// Calculate the price of the European option
double EuropeanOption::Price() const
{
	double d1 = D1();
	double d2 = D2(d1);

	if (m_type == "Call")
	{
		return m_S * std::exp((m_b - m_r) * m_T) * N(d1) - m_K * std::exp(-m_r * m_T) * N(d2);
	}

	if (m_type == "Put")
	{
		return m_K * std::exp(-m_r * m_T) * (1 - N(d2)) - m_S * std::exp((m_b - m_r) * m_T) * (1 - N(d1));
	}
}

// Calculate the price using put-call parity
double EuropeanOption::PricePutCallParity() const
{
	if (m_type == "Call")
	{
		return Price() + m_K * std::exp(-m_r * m_T) - m_S;
	}

	if (m_type == "Put")
	{
		return Price() + m_S - m_K * std::exp(-m_r * m_T);
	}
}

// Calculate the Delta of the option
double EuropeanOption::Delta() const
{
	if (m_type == "Call")
	{
		return exp((m_b - m_r) * m_T) * N(D1());
	}

	if (m_type == "Put")
	{
		return -exp((m_b - m_r) * m_T) * (1 - N(D1()));
	}
}

// Calculate the Gamma of the option
double EuropeanOption::Gamma() const
{
	return N_prime(D1()) * exp((m_b - m_r) * m_T) / (m_S * m_sigma * sqrt(m_T));
}

// Calculate the Vega of the option
double EuropeanOption::Vega() const
{
	return m_S * sqrt(m_T) * exp((m_b - m_r) * m_T) * N_prime(D1());
}

// Calculate the Theta of the option
double EuropeanOption::Theta() const
{
	double d1 = D1();
	double d2 = D2(d1);

	if (m_type == "Call")
	{
		return (-(m_S)*m_sigma * exp((m_b - m_r) * m_T) * N_prime(d1) / (2 * sqrt(m_T))) - ((m_b - m_r) * m_S * exp((m_b - m_r) * m_T) * N(d1)) - (m_r * m_K * exp(-m_r * m_T) * N(d2));
	}

	if (m_type == "Put")
	{
		return (-(m_S)*m_sigma * exp((m_b - m_r) * m_T) * N_prime(d1) / (2 * sqrt(m_T))) - ((m_b - m_r) * m_S * exp((m_b - m_r) * m_T) * (1 - N(d1))) - (m_r * m_K * exp(-m_r * m_T) * (1 - N(d2)));
	}
}

// Calculate the Rho of the option
double EuropeanOption::Rho() const
{
	double d1 = D1();
	double d2 = D2(d1);

	if (m_type == "Call")
	{
		return m_K * m_T * exp(-m_r * m_T) * N(d2);
	}

	if (m_type == "Put")
	{
		return -m_K * m_T * exp(-m_r * m_T) * (1 - N(d2));
	}
}

// Calculate the Vanna of the option
double EuropeanOption::Vanna() const
{
	double d1 = D1();

	return exp((m_b - m_r) * m_T) * N_prime(d1) * (d1 / m_sigma - 1) * m_S * sqrt(m_T);
}

// Calculate the Charm of the option
double EuropeanOption::Charm() const
{
	double d1 = D1();
	double d2 = D2(d1);
	double term1 = -exp((m_b - m_r) * m_T) * N_prime(d1) * (2 * (m_r - m_b) * m_T - d2 * m_sigma * sqrt(m_T)) / (2 * m_T * m_sigma * sqrt(m_T));

	if (m_type == "Call")
	{
		return term1 - (m_r - m_b) * exp((m_b - m_r) * m_T) * N(d1);
	}
	if (m_type == "Put")
	{
		return term1 + (m_r - m_b) * exp((m_b - m_r) * m_T) * (1 - N(d1));
	}
}

// Calculate the Speed of the option
double EuropeanOption::Speed() const
{
	double d1 = D1();

	return -N_prime(d1) * exp((m_b - m_r) * m_T) * (d1 / (m_S * m_S * m_sigma * sqrt(m_T)));
}

// Calculate the Color of the option
double EuropeanOption::Color() const
{
	double d1 = D1();

	return -N_prime(d1) * exp((m_b - m_r) * m_T) * (2 * (m_r - m_b) * m_T - d1 * m_sigma * sqrt(m_T)) / (2 * m_T * m_S * m_sigma * sqrt(m_T));
}

// Calculate the DvegaDtime of the option
double EuropeanOption::DvegaDtime() const
{
	double d1 = D1();

	return m_S * sqrt(m_T) * exp((m_b - m_r) * m_T) * N_prime(d1) * (m_r - m_b - d1 * m_sigma / (2 * sqrt(m_T)));
}

// Calculate the Vomma of the option
double EuropeanOption::Vomma() const
{
	double d1 = D1();

	return Vega() * d1 * (d1 - 1) / m_sigma;
}

// Calculate the Veta of the option
double EuropeanOption::Veta() const
{
	double d1 = D1();

	return -m_S * exp((m_b - m_r) * m_T) * N_prime(d1) * sqrt(m_T) * (m_r - m_b + (d1 * m_sigma) / (2 * sqrt(m_T)));
}

// Calculate the Zomma of the option
double EuropeanOption::Zomma() const
{
	double d1 = D1();

	return Gamma() * (d1 * (d1 - 1) - 1) / m_sigma;
}

// Calculate the Lambda of the option
double EuropeanOption::Lambda() const
{
	return Delta() * (m_S / Price());
}

// Calculate the Ultima of the option
double EuropeanOption::Ultima() const
{
	double d1 = D1();
	double d2 = D2(d1);

	return -Vega() * (d1 * (d1 - 3) * (d1 - 1) - 1) / (m_sigma * m_sigma);
}

// Calculate the numeric Delta of the option
double EuropeanOption::NumericDelta(const double& h) const
{
	double up_price = EuropeanOption(m_type, m_T, m_K, m_S + h, m_r, m_sigma, m_id, m_b).Price();

	double down_price = EuropeanOption(m_type, m_T, m_K, m_S - h, m_r, m_sigma, m_id, m_b).Price();

	return (up_price - down_price) / (2 * h);
}

// Calculate the numeric Gamma of the option
double EuropeanOption::NumericGamma(const double& h) const
{
	double up_price = EuropeanOption(m_type, m_T, m_K, m_S + h, m_r, m_sigma, m_id, m_b).Price();

	double down_price = EuropeanOption(m_type, m_T, m_K, m_S - h, m_r, m_sigma, m_id, m_b).Price();

	return (up_price - 2 * Price() + down_price) / (h * h);
}

// Check if put-call parity is satisfied
void EuropeanOption::CheckPutCallParity(const double& market_price, const double& threshold) const
{
	double implied_price = PricePutCallParity();

	double rel_diff = std::abs(market_price - implied_price) / implied_price;

	if (rel_diff <= threshold)
	{
		std::cout << "Put-Call Parity is satisfied within a " << threshold << " threshold" << std::endl;
	}
	else
	{
		std::cout << "Put-Call Parity is not satisfied within a " << threshold << " threshold" << std::endl;
	}
}

// Set the option type
EuropeanOption& EuropeanOption::type(const std::string& type)
{
	m_type = type;
	return *this;
}

// Set the option maturity time
EuropeanOption& EuropeanOption::T(const double& T)
{
	m_T = T;
	return *this;
}

// Set the option strike price
EuropeanOption& EuropeanOption::K(const double& K)
{
	m_K = K;
	return *this;
}

// Set the option spot price
EuropeanOption& EuropeanOption::S(const double& S)
{
	m_S = S;
	return *this;
}

// Set the risk-free interest rate
EuropeanOption& EuropeanOption::r(const double& r)
{
	m_r = r;
	return *this;
}

// Set the volatility of the option
EuropeanOption& EuropeanOption::sigma(const double& sigma)
{
	m_sigma = sigma;
	return *this;
}

// Set the cost of carry
EuropeanOption& EuropeanOption::b(const double& b)
{
	m_b = b;
	return *this;
}

// Define << ostream operator function
std::ostream& operator << (std::ostream& os, const EuropeanOption& source)
{
	// Sends description to output stream
	os << "Option " << source.m_id << ": " << source.m_type << ", T: " << source.m_T << ", K: " << source.m_K << ", S: "
		<< source.m_S << ", r: " << source.m_r << ", sigma: " << source.m_sigma << ", b: " << source.m_b;

	// Returns the output stream
	return os;
}
