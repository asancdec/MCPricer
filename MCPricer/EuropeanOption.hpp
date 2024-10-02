// (C++) Monte Carlo Option Pricer with Euler - Maruyama Discretization
// EuropeanOption.hpp
// Álvaro Sánchez de Carlos
// Description: this file contains the header code for the EuropeanOption class

// If EUROPEANOPTION_HPP is not defined
#ifndef EUROPEANOPTION_HPP
// Define EUROPEANOPTION_HPP
#define EUROPEANOPTION_HPP

#include <string>
#include <ostream>
#include <vector>

// Class definition for EuropeanOption
class EuropeanOption
{
private:

	// Option type (call or put)
	std::string m_type;
	// Time to expiration
	double m_T;
	// Strike price
	double m_K;
	// Spot price
	double m_S;
	// Risk-free interest rate
	double m_r;
	// Volatility
	double m_sigma;
	// Cost of carry
	double m_b;
	// Option ID
	int m_id;

	// Calculate d1 for Black-Scholes formula
	double D1() const;
	// Calculate d2 for Black-Scholes formula
	double D2(const double& d1) const;

	// Standard normal cumulative distribution function
	double N(double x) const;
	// Standard normal probability density function
	double N_prime(double x) const;

public:

	// Constructor
	EuropeanOption(const std::string& type, const double& T, const double& K, const double& S, const double& r,
		const double& sigma, const int& id = 1, const double& b = std::numeric_limits<double>::quiet_NaN());
	// Copy constructor
	EuropeanOption(const EuropeanOption& source);
	// Assignment operator
	EuropeanOption& operator=(const EuropeanOption& source);
	// Convert option parameters to vector of strings
	std::vector<std::string> ConvertToVectorString() const;

	// Pricing functions
	// Calculate option price
	double Price() const;
	// Calculate price using Put-Call Parity
	double PricePutCallParity() const;

	// Greeks functions
	// Calculate Delta
	double Delta() const;
	// Calculate Gamma
	double Gamma() const;
	// Calculate Vega
	double Vega() const;
	// Calculate Theta
	double Theta() const;
	// Calculate Rho
	double Rho() const;
	// Calculate Vanna
	double Vanna() const;
	// Calculate Charm
	double Charm() const;
	// Calculate Speed
	double Speed() const;
	// Calculate Color
	double Color() const;
	// Calculate DvegaDtime
	double DvegaDtime() const;
	// Calculate Vomma
	double Vomma() const;
	// Calculate Veta
	double Veta() const;
	// Calculate Zomma
	double Zomma() const;
	// Calculate Lambda
	double Lambda() const;
	// Calculate Ultima
	double Ultima() const;

	// Numeric functions
	// Calculate Numeric Delta
	double NumericDelta(const double& h) const;
	// Calculate Numeric Gamma
	double NumericGamma(const double& h) const;

	// Check Put-Call Parity
	void CheckPutCallParity(const double& market_price, const double& threshold = 0.05) const;

	// Parameter modification functions
	// Set option type
	EuropeanOption& type(const std::string& type);
	// Set time to expiration
	EuropeanOption& T(const double& T);
	// Set strike price
	EuropeanOption& K(const double& K);
	// Set spot price
	EuropeanOption& S(const double& S);
	// Set risk-free interest rate
	EuropeanOption& r(const double& r);
	// Set volatility
	EuropeanOption& sigma(const double& sigma);
	// Set cost of carry
	EuropeanOption& b(const double& b);

	// Get inline functions
	// Get option type
	const std::string& type() const { return m_type; }
	// Get time to expiration
	const double& T() const { return m_T; }
	// Get strike price
	const double& K() const { return m_K; }
	// Get spot price
	const double& S() const { return m_S; }
	// Get risk-free interest rate
	const double& r() const { return m_r; }
	// Get volatility
	const double& sigma() const { return m_sigma; }
	// Get cost of carry
	const double& b() const { return m_b; }
	// Get option ID
	const int& id() const { return m_id; }

	// Friend functions
	// Define << ostream operator function
	friend std::ostream& operator << (std::ostream& os, const EuropeanOption& source);

};

// End of the conditional inclusion of the header file
#endif
