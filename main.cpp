/*
 * main.cpp
 *
 *  Created on: Oct 6, 2017
 *      Author: ricardo
 */

#include <fstream>
#include "Utils.hpp"
#include "EvoStrat.hpp"

int main() {

    inputs in;
    in = readParameters("input_parameters.txt");
    in.print();

    //initializatrion
    std::vector<double> _range;
    _range.push_back(in.max);
    _range.push_back(in.min);

    if (in.tau == -1) {
        in.tau = 1/std::sqrt(in.n_var);
    }

    if (in.tau_prime == -1) {
        in.tau_prime = 1/std::sqrt(2* std::sqrt(in.n_var));
    }

    std::vector<individual> best;

    population lambda_population;
    lambda_population.mu = in.lambda;

    std::vector<double> best_sol;
    std::vector<int> generation;
    std::vector<std::vector<double>> sims;

    double (*test_function)(std::vector<double>);

    if (in.function == 1) {
        test_function = sphere;
    }

    if (in.function == 2) {
        test_function = schwefel;
    }

    if (in.function == 3) {
        test_function = deJong;
    }

    //test parameters
    std::vector<double> test_sigmas {0.002,0.02,0.2,2.0,20.0,200.0};

    std::vector<double> test_taus {0.01/std::sqrt(in.n_var),
    	0.1/std::sqrt(in.n_var),
		0.5/std::sqrt(in.n_var),
		1/std::sqrt(in.n_var),
		2/std::sqrt(in.n_var),
		3/std::sqrt(in.n_var)};

    std::vector<double> test_tau_prime {std::sqrt(0.5*test_taus[0]),
    	std::sqrt(0.5*test_taus[1]),
		std::sqrt(0.5*test_taus[2]),
		std::sqrt(0.5*test_taus[3]),
		std::sqrt(0.5*test_taus[4]),
		std::sqrt(0.5*test_taus[5])};

    std::vector<double> test_epsilon {0.001,0.01,0.1,1.0,10.0,100.0};

    if (in.test_parameters != 0) {

        //adapt simulations to test size
        if (in.test_parameters == 1 || in.test_parameters == 3) {
            in.simulations = test_sigmas.size();
        } else {
            in.simulations = test_taus.size();
        }
    }

    int test_n = 0;  //initialize test number

    //
    //evolution strategy
    for (int nsims=0; nsims < in.simulations; nsims++){

        if (in.test_parameters != 0 ) {
            if (in.test_parameters == 1) { //test sigma
                in.sigma = test_sigmas[test_n];
                test_n++;
            }

            if (in.test_parameters == 2) { //test tau
                in.tau = test_taus[test_n];
                in.tau_prime = test_tau_prime[test_n];
                test_n++;
            }

            if (in.test_parameters == 3) { //test epsilon
                in.epsilon = test_epsilon[test_n];
                test_n++;
            }
        }

        //initialize
        population pob = population_init(test_function, _range,
        		in.n_var, in.mu, in.sigma);

        for (int t=0; t < in.generations; t++){
            //generate lambda childs
            for (int i = 0; i < in.lambda; i++){

                individual lambda_child = combination(test_function,
                		pob, in.rho, in.type_combination);

                mutation(test_function, lambda_child, in.tau,
                		in.tau_prime, in.type_mutation, in.epsilon, _range);

                lambda_population.individuals.push_back(lambda_child);

                lambda_child.var.clear();
                lambda_child.sigma.clear();

            } //for lambda

            //population is return already sorted
            selection(pob,lambda_population,in.type_selection);

            lambda_population.individuals.clear();

            // store values
            best.push_back(pob.individuals[0]);
            best_sol.push_back(best.back().fitness);
            generation.push_back(t);

        }//for generations

        sims.push_back(best_sol);

        if (in.test_parameters == 0) {
            std::cout<<"*** End sim " << nsims + 1 << " ***\n";
        } else {
            std::cout<<"*** End test " << nsims + 1 << " ***\n";
            std::cout <<"Sigma     : " << in.sigma << "\n";
            std::cout <<"Tau       : " << in.tau << "\n";
            std::cout <<"Tau_prime : " << in.tau_prime << "\n";
            std::cout <<"Epsilon   : " << in.epsilon << "\n\n";
        }
        //reset
        best_sol.clear();
        generation.clear();
        pob.individuals.clear();
        lambda_population.individuals.clear();

    }// simulations

    //output
    std::ofstream F;
    F.open("output.csv");

    for (int i = 0; i < in.generations; i ++) {
        F << generation[i] <<",";
        for (int j=0;j < in.simulations; j++) {
            F << sims[j][i] <<",";
        }
        F <<"\n";
    }
    F.close();
    std::cout<<"*** FIN: Simulations can be found in output.csv ***\n";
    return 0;

}





