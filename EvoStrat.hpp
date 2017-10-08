/*
 * EvoStrat.hpp
 *
 *  Created on: Oct 8, 2017
 *      Author: ricardo
 */

#ifndef EVOSTRAT_HPP_
#define EVOSTRAT_HPP_

#include <random>
#include <vector>
#include <iostream>
#include <algorithm>

struct individual {
    std::vector<double> var;
    std::vector<double> sigma;
    double fitness;
};

struct population {
    std::vector<individual> individuals;
    int mu;

};

void mutation(double (*fp)(std::vector<double>),
		individual& x,
		double tau,
		double tau_prime,
		int type,
		double precision,
		std::vector<double> _range);

population population_init(double (*fp)(std::vector<double>),
		std::vector<double> _range,
		int n_var,
		int mu,
		double sigma);

individual combination(double (*fp)(std::vector<double>),
		population p,
		int rho,
		int type);

void selection(population& parent_population,
		population lambda_population,
		int type);

template<class T> void print_vector(std::vector<T> v) {
    std::cout <<"\n[";
    for(auto it = v.begin(); it != v.end();++it){
        std::cout << (*it) << " , ";
    }
    std::cout << "]\n";

}





#endif /* EVOSTRAT_HPP_ */
