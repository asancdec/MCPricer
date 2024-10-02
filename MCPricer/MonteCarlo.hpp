// (C++) Monte Carlo Option Pricer with Euler - Maruyama Discretization
// MonteCarlo.hpp
// Álvaro Sánchez de Carlos
// Description: this file contains the header code of the derived MonteCarlo class

// If MONTECARLO_HPP is not defined
#ifndef MONTECARLO_HPP
// Define MONTECARLO_HPP
#define MONTECARLO_HPP

#include "EuropeanOption.hpp"

// Define MonteCarlo derived class from EuropeanOption
class MonteCarlo : public EuropeanOption
{
private:

    // Declare private member variables
    long m_subintervals;
    long m_simulations;

    // Declare SD private function
    double SD(const double& sum_payoff, const double& sum_square_payoff) const;
    // Declare SE private function
    double SE(const double& sd) const;

public:

    // Constructor 
    MonteCarlo(const EuropeanOption& option, const long& subintervals = 1e2, const long& simulations = 1e4);

    // Copy constructor
    MonteCarlo(const MonteCarlo& source);
    
    // Assignement operator
    MonteCarlo& operator=(const MonteCarlo& source);

    // Declare the Price function
    double Price(const double& beta = 1, const bool& error_analysis = true) const;
};

// End of the conditional inclusion of the header file
#endif