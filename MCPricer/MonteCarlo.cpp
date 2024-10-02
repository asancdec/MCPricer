//(C++) Monte Carlo Option Pricer with Euler - Maruyama Discretization
// MonteCarlo.cpp
// Álvaro Sánchez de Carlos
// Description: this file contains the source code of the derived MonteCarlo class

#include <algorithm>
#include <functional>
#include <iomanip>
#include <boost/random.hpp>
#include "MonteCarlo.hpp"

// Define SD private function
double MonteCarlo::SD(const double& sum_payoff, const double& sum_square_payoff) const
{
    return std::sqrt(((sum_square_payoff - (sum_payoff * sum_payoff / m_simulations)) / (m_simulations - 1))) * exp(this->r() * this->T());
}

// Define SE private function
double MonteCarlo::SE(const double& sd) const
{
    return sd / std::sqrt(m_simulations);
}

// Constructor 
MonteCarlo::MonteCarlo(const EuropeanOption& option, const long& subintervals, const long& simulations) :
    EuropeanOption(option),
    m_subintervals(subintervals),
    m_simulations(simulations)
{}

// Copy Constructor
MonteCarlo::MonteCarlo(const MonteCarlo& source) :
    EuropeanOption(source),
    m_subintervals(source.m_subintervals),
    m_simulations(source.m_simulations)
{}

// Assignment operator
MonteCarlo& MonteCarlo::operator=(const MonteCarlo& source)
{
    // Check for self assignment
    if (this == &source)
        return *this;

    EuropeanOption::operator=(source);
    m_subintervals = source.m_subintervals;
    m_simulations = source.m_simulations;

    return *this;
}

// Define the Price function
double MonteCarlo::Price(const double& beta, const bool& error_analysis) const
{   
    // Extract option parameters
    const double K = this->K();
    const double T = this->T();
    const double r = this->r();
    const double sigma = this->sigma();
    const double S = this->S();
    const std::string type = this->type();

    // Precompute MC parameters to optimize speed
    const double tn = T / m_subintervals;
    const double drift_const = r * tn;
    const double diffusion_const = sigma * std::sqrt(tn);

    // Create a random number generator (Mersenne Twister)
    boost::random::mt19937 wiener_process;

    // Define a normal distribution
    boost::random::normal_distribution<> dist(0, 1);

    // Define a variable to store payoffs
    double sum_payoff = 0.0;

    // Define a variable to store squared payoffs
    double sum_square_payoff = 0.0;

    // Declare double SN variable
    double SN;

    // Precompute beta-specific logic outside loops to avoid runtime branching
    auto pow_function = [beta](double x) 
        {
        if (beta == 1.0) return x;
        else if (beta == 0.5) return std::sqrt(x);
        else if (beta == 2.0) return x * x;
        else return std::pow(x, beta);
        };

    // Parallelize the loop over m_simulations
       #pragma omp parallel reduction(+:sum_payoff, sum_square_payoff)
    {
        // Thread-local random number generator
        boost::random::mt19937 local_wiener_process(wiener_process);
        #pragma omp for
        for (long i = 1; i <= m_simulations; ++i)
        {
            // Re start the S0 variable to the current underlying spot price (S) every simulation
            double S0 = S;

            // Simulate one path in the underlying for N subintervals
            for (long a = 1; a <= m_subintervals; ++a)
            {
                // Avoid unnecessary std::pow calls
                SN = S0 + drift_const * S0 + diffusion_const * pow_function(S0) * dist(local_wiener_process);

                // Update the initial Spot Value for the next subinterval
                S0 = SN;
            }

            // Calculate the payoff at the end of the simulation
            double payoff = (type == "Call")
                ? std::max(SN - K, 0.0)
                : std::max(K - SN, 0.0);

            // Add the payoff to accumulated payoffs
            sum_payoff += payoff;

            // Add the squared payoff to accumulated squared payoffs
            if (error_analysis) sum_square_payoff += (payoff * payoff);
        }
    }

    // Calculate the average payoff
    const double average_payoff = sum_payoff / m_simulations;

    // Calculate the price with exponential discount of the average payoff
    double price = average_payoff * exp(-r * T);

    // If error_analysis, print an analysis of the errors
    if (error_analysis)
    {   
        // Calculate SD
        double sd = SD(sum_payoff, sum_square_payoff);

        // Print the error results as a table
        std::cout << std::setw(20) << "Simulations"
            << std::setw(20) << "Subintervals"
            << std::setw(20) << "BSM Price"
            << std::setw(20) << "MC Price"
            << std::setw(20) << "SD"
            << std::setw(20) << "SE"
            << std::endl;

        std::cout << std::setw(20) << m_simulations
            << std::setw(20) << m_subintervals
            << std::setw(20) << this->EuropeanOption::Price()
            << std::setw(20) << price
            << std::setw(20) << sd
            << std::setw(20) << SE(sd)
            << std::endl;
    }

    // Return the price
    return price;
}
