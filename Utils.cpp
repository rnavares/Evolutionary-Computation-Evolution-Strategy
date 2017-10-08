/*
 * Utils.cpp
 *
 *  Created on: Oct 8, 2017
 *      Author: ricardo
 */
#include "Utils.hpp"


double sphere(std::vector<double> vars) {

    double res = 0;
    for (auto it = vars.begin(); it != vars.end(); ++it){
       res = res + (*it) * (*it);
    }
    return res;
}

double schwefel(std::vector<double> vars) {

    double res = 0;
    for (auto it = vars.begin(); it != vars.end(); ++it){
       //res = res - std::sin(*it) * std::sin((*it)*(*it));
       res = res - (*it) * std::sin(std::sqrt(std::abs((*it))));
    }

    res = res + 418.9829 * vars.size();

    return res;
}

double deJong(std::vector<double> vars) {


    std::vector<std::vector<double>> A;

    std::vector<double> row {-32,-16,0,16,32};
    std::vector<double> A1;
    for(int i = 0; i < 5; i++) {
        A1.insert(A1.end(),row.begin(),row.end());
    }

    std::for_each(row.begin(),row.end(), [](double &x){x=-32;});

    std::vector<double> A2;
    for(int i = 0; i < 5; i++) {
        A2.insert(A2.end(),row.begin(),row.end());
        std::for_each(row.begin(),row.end(), [](double &x){x=x+16;});
    }

    A.push_back(A1);
    A.push_back(A2);

    double res = 0.0;

    double sum_a = 0.0;
    double sum_x = 0.0;
    for (int i = 0; i < A1.size(); i++) {
        for(int j = 0; j < vars.size(); j++) {
            sum_x = sum_x + std::pow(vars[j] - A[j][i],6.0);
        }
        sum_a = sum_a + 1/((i+1) + sum_x);
        sum_x = 0;
    }

    sum_a = sum_a + 1/500;

    res = 1/sum_a;

    return res;
}

inputs readParameters(std::string F) {
    std::ifstream fid;

    fid.open(F);

    inputs e;

    while(!fid.eof()){
        std::string dat;
        std::string symbol;
        std::string parameter;



        fid >> parameter;

        if(parameter.compare("max") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.max;
        }

        if(parameter.compare("min") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.min;
        }

        if(parameter.compare("n_var") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.n_var;
        }

        if(parameter.compare("mu") == 0) {

            fid >> symbol >> dat;
            std::istringstream(dat) >> e.mu;
        }

        if(parameter.compare("lambda") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.lambda;
        }

        if(parameter.compare("sigma") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.sigma;
        }

        if(parameter.compare("rho") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.rho;
        }

        if(parameter.compare("epsilon") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.epsilon;
        }

        if(parameter.compare("tau") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.tau;
        }

        if(parameter.compare("tau_prime") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.tau_prime;
        }

        if(parameter.compare("type_mutation") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.type_mutation;
        }

        if(parameter.compare("type_combination") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.type_combination;
        }

        if(parameter.compare("type_selection") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.type_selection;
        }

        if(parameter.compare("generations") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.generations;
        }

        if(parameter.compare("simulations") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.simulations;
        }

        if(parameter.compare("function") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.function;
        }

        if(parameter.compare("test_parameters") == 0) {
            fid >> symbol >> dat;
            std::istringstream(dat) >> e.test_parameters;
        }
    }

    fid.close();

    return(e);
}


