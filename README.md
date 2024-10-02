# Monte Carlo Option Pricer

## Description

The **Monte Carlo Option Pricer** is a C++ implementation of a financial option pricing model using the Monte Carlo simulation method with the Euler-Maruyama discretization technique. The program is designed to price European call and put options by simulating random paths of the underlying asset price and calculating the expected payoff at maturity. 

The model uses the Black-Scholes framework to derive parameters for the options and performs error analysis to compare the Monte Carlo price against the theoretical Black-Scholes price.

## Features

- **Monte Carlo Simulation**: Implements the Monte Carlo method to estimate the option prices.
- **Euler-Maruyama Discretization**: Uses Euler-Maruyama for simulating paths of the underlying asset.
- **Error Analysis**: Optional error analysis providing standard deviation and standard error of the estimated prices.
- **European Options**: Specifically designed for European-style options (call and put).
- **Boost Library Integration**: Utilizes the Boost library for random number generation and statistical distributions.

## Files

- `EuropeanOption.hpp`: Header file containing the declaration of the `EuropeanOption` class, which defines the properties and pricing methods for European options.
- `EuropeanOption.cpp`: Contains the implementation of the `EuropeanOption` class.
- `MonteCarlo.hpp`: Header file containing the declaration of the `MonteCarlo` class, which performs the Monte Carlo simulation to price options.
- `MonteCarlo.cpp`: Contains the implementation of the `MonteCarlo` class.
- `MCPricer.cpp`: The main driver program that creates instances of `EuropeanOption` and `MonteCarlo`, runs simulations, and displays results.

## Usage

1. **Compile the Code**: Use a C++ compiler (e.g., g++) to compile the source files. Make sure to link against the Boost library. 

   ```bash
   g++ -o MonteCarloOptionPricer MCPricer.cpp EuropeanOption.cpp MonteCarlo.cpp -lboost_system -lboost_random
