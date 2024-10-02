// (C++) Monte Carlo Option Pricer with Euler - Maruyama Discretization
// MCPricer.cpp
// Álvaro Sánchez de Carlos
// Description: This file contains the main function of the MCPricer

#include <iostream>
#include <vector>
#include "EuropeanOption.hpp"
#include "MonteCarlo.hpp"

// Define main function of the program
int main()
{
    // Create a European call option with specified parameters
    EuropeanOption call_option("Call", 0.25, 65, 60, 0.08, 0.3, 1);
    // Print the details of the call option
    std::cout << call_option << std::endl;

    // Initialize the number of simulations
    int n_simulations = 10;
    // Loop to increase the number of simulations and price the option
    for (int i = 1; i <= 7; ++i)
    {
        // Create a Monte Carlo simulation with the call option, 100 subintervals, and n_simulations
        // Price the option and print the results
        MonteCarlo(call_option, 100, n_simulations).Price(1, true);
        // Increase the number of simulations by a factor of 10
        n_simulations *= 10;
    }

    // Initialize the number of subintervals
    int n_subintervals = 10;
    // Loop to increase the number of subintervals and price the option
    for (int i = 1; i <= 7; ++i)
    {
        // Create a Monte Carlo simulation with the call option, n_subintervals subintervals, and 1000 simulations
        // Price the option and print the results
        MonteCarlo(call_option, n_subintervals, 1000).Price(1, true);
        // Increase the number of subintervals by a factor of 10
        n_subintervals *= 10;
    }

    // Create a European put option with specified parameters
    EuropeanOption put_option("Put", 1.0, 100, 100, 0.00, 0.2, 2);
    // Print the details of the put option
    std::cout << put_option << std::endl;

    // Initialize the number of simulations
    n_simulations = 10; // Reset for put option
    // Loop to increase the number of simulations and price the option
    for (int i = 1; i <= 7; ++i)
    {
        // Create a Monte Carlo simulation with the put option, 100 subintervals, and n_simulations
        // Price the option and print the results
        MonteCarlo(put_option, 100, n_simulations).Price(1, true);
        // Increase the number of simulations by a factor of 10
        n_simulations *= 10;
    }

    // Initialize the number of subintervals
    n_subintervals = 10; // Reset for put option
    // Loop to increase the number of subintervals and price the option
    for (int i = 1; i <= 7; ++i)
    {
        // Create a Monte Carlo simulation with the put option, n_subintervals subintervals, and 1000 simulations
        // Price the option and print the results
        MonteCarlo(put_option, n_subintervals, 1000).Price(1, true);
        // Increase the number of subintervals by a factor of 10
        n_subintervals *= 10;
    }

    // Return 0 to indicate successful execution
    return 0;
}
