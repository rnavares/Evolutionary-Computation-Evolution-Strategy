/*
 * Utils.hpp
 *
 *  Created on: Oct 8, 2017
 *      Author: ricardo
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>

struct inputs {
    double max = 0;              //search space
    double min = 0;
    int mu = 0;                  //initial population size
    int n_var = 0;               //number of variables
    int type_mutation = 0;       //1: Uncorr 1 Step  2: Uncorr. N step
    int type_combination = 0;    //1: Discreta  2: intermedia
    double sigma = 0.0;
    int lambda = 0;
    int rho = 0;                 //number of parents
    double epsilon = 0.0;
    double tau = 0;              // learning rate, -1: default
    double tau_prime = 0;
    int type_selection = 0;      //1: (mu,lambda)    2: (mu+lambda)
    int generations = 0;         //stop condition
    int simulations = 0;         //time the algorithm will be executed
    int function = 0;           // 1: Sphere   2: Schwefel  3: DeJong
    int test_parameters = 0;

    void print() {
        std::cout << "\nPARAMETROS";
        std::cout << "\n Rango [" << max <<","<< min <<"]";
        std::cout << "\n Numero de variables    :\t" << n_var;
        std::cout << "\n Mu                     :\t" << mu;
        std::cout << "\n Lambda                 :\t" << lambda;
        std::cout << "\n Sigma                  :\t"<< sigma;
        std::cout << "\n Rho                    :\t" << rho;
        std::cout << "\n Epsilon                :\t"<< epsilon;
        std::cout << "\n Tau                    :\t" << tau;
        std::cout << "\n Tau prima              :\t" << tau_prime;
        std::cout << "\n Tipo Mutacion          :\t" << type_mutation;
        std::cout << "\n Tipo de Recombinacion  :\t" << type_combination;
        std::cout << "\n Tipo de Seleccion      :\t" << type_selection;
        std::cout << "\n Generaciones           :\t" << generations;
        std::cout << "\n Funcion de Test        :\t" << function;
        std::cout << "\n Simulaciones           :\t" << simulations;
        std::cout << "\n\n";

    }

};


//funcion Sphere
double sphere(std::vector<double> vars);
//funcion de Schwefel
double schwefel(std::vector<double> vars);
//funcion de DeJong5
double deJong(std::vector<double> vars);
inputs readParameters(std::string F);




#endif /* UTILS_HPP_ */
