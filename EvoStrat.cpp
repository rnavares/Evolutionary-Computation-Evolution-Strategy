/*
 * EvoStrat.cpp
 *
 *  Created on: Oct 8, 2017
 *      Author: ricardo
 */

#include "EvoStrat.hpp"

population population_init(double (*fp)(std::vector<double>),
		std::vector<double> _range,
		int n_var,
		int mu,
		double sigma) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(_range[1],_range[0]);

    individual xi;
    population pob;

    pob.mu = mu;

    for (int i = 0 ; i < mu; i++) {
        for (int j = 0 ; j < n_var; j++) {
            xi.var.push_back(dist(gen));
            xi.sigma.push_back(sigma);
        }
        xi.fitness = fp(xi.var); //fitness
        pob.individuals.push_back(xi);
        xi.var.clear();
        xi.sigma.clear();

    }

    return pob;
}


individual combination(double (*fp)(std::vector<double>),
		population p,
		int rho,
		int type) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0,p.mu-1);
    std::uniform_int_distribution<> dist2(0,rho-1);

    individual child;

    //generate rho candidates
    std::vector<int> parents;
    for (int i = 0; i < rho; i++) {
        parents.push_back(dist(gen));   //no controlamos mismo padre
    }

    int n_variables = p.individuals[0].var.size();

    if (type == 1) {  //discrete combination

        for(int i = 0; i < n_variables; i++){   // for each variable we select random candidates
                                                // among the parents to transfer to the children
            int _index = parents[dist2(gen)];
            child.var.push_back(p.individuals[_index].var[i]);
            child.sigma.push_back(p.individuals[_index].sigma[i]);
        }
        child.fitness = fp(child.var);

    }

    if (type == 2) {  // intermediate combination

        for(int i = 0; i < n_variables; i++){   // transfer to the child the average
                                                // of the values of the parents
            double sum_xi = 0.0;
            double sum_VARi = 0.0;

            for (int j = 0; j < rho; j++) {
                sum_xi = sum_xi + p.individuals[parents[j]].var[i] ;
                sum_VARi = sum_VARi + p.individuals[parents[j]].sigma[i] *
                		p.individuals[parents[j]].sigma[i]; //aggregate variaces
            }
            child.var.push_back(sum_xi / rho);
            child.sigma.push_back(std::sqrt(sum_VARi / rho));
        }
        child.fitness = fp(child.var);

    }

   // print_vector(child.var);
   // print_vector(child.sigma);
   // std::cout << "FIT " << child.fitness;


    return child;
}


void mutation(double (*fp)(std::vector<double>),
		individual& x,
		double tau,
		double tau_prime,
		int type,
		double precision,
		std::vector<double> _range){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0,1.0);

    if (type == 1) { // Uncorr. w 1 Step

       // generate disturbance, only one for all sigmas
       double disturbance = std::exp(tau*dist(gen));

       for (int i = 0; i < (int) x.var.size(); i++) {

            x.sigma[i] = x.sigma[i] * disturbance;

            if (x.sigma[i] < precision) {
                x.sigma[i] = precision;
            }

            x.var[i] = x.var[i] + x.sigma[i] * dist(gen);

            //controlamos los limites
            if (x.var[i] > _range[0]) {
                x.var[i] = _range[0];
            }

            if (x.var[i] < _range[1]) {
                x.var[i] = _range[1];
            }

       }

       x.fitness = fp(x.var);


    }

    if (type == 2) { // Uncorr. w N Step

        double N01 = dist(gen);
        for (int i = 0; i < (int) x.var.size(); i++) {
            //add disturbance to each sigma
            double disturbance = std::exp(tau_prime * N01 + tau * dist(gen));

            x.sigma[i] = x.sigma[i] * disturbance;

            if (x.sigma[i] < precision) {
                x.sigma[i] = precision;
            }

            dist.reset();

            x.var[i] = x.var[i] + x.sigma[i] * dist(gen);

            //control los limits
            if (x.var[i] > _range[0]) {
                x.var[i] = _range[0];
            }

            if (x.var[i] < _range[1]) {
                x.var[i] = _range[1];
            }

       }

        x.fitness = fp(x.var);
    }

}

void selection(population& parent_population,
		population lambda_population,
		int type) {

    if (type == 1) { // (mu,lambda)

        if (parent_population.mu >= lambda_population.mu) {
            std::cout <<"\n ERROR: mu has to be < lambda \n";
            throw(1);
        }

        std::vector<std::vector<double>> fitness;
        //enumerate a vector with fit and _index
        for (int i=0; i < lambda_population.mu; i++) {
            std::vector<double> row;
            row.push_back(i);
            row.push_back(lambda_population.individuals[i].fitness);
            fitness.push_back(row);

        }


        //sort by fitness
        sort(fitness.begin(),fitness.end(),
        		[](std::vector<double> f1,
        				std::vector<double>f2)->bool{return (f1[1]<f2[1]);});


        //replace best lambdas in the initial population
        for (int i = 0; i < parent_population.mu;i++){
            parent_population.individuals[i] = lambda_population.individuals[fitness[i][0]];
        }

        fitness.clear();
    }

    if (type == 2) { //(mu+lambda)

        std::vector<std::vector<double>> fitness_lambda;
        //enumerate a vector with fit and _index
        for (int i=0; i < lambda_population.mu; i++) {
            std::vector<double> row;
            row.push_back(i);
            row.push_back(lambda_population.individuals[i].fitness);
            fitness_lambda.push_back(row);
        }

        //sort by fitness
        sort(fitness_lambda.begin(),fitness_lambda.end(),
        		[](std::vector<double> f1,
        				std::vector<double>f2)->bool{return (f1[1]<f2[1]);});

        std::vector<std::vector<double>> fitness_mu;
        //enumerate a vector with fit and _index
        for (int i=0; i < parent_population.mu; i++) {
            std::vector<double> row;
            row.push_back(i);
            row.push_back(parent_population.individuals[i].fitness);
            fitness_mu.push_back(row);
        }

        //sort by fitness
        sort(fitness_mu.begin(),fitness_mu.end(),
        		[](std::vector<double> f1,
        				std::vector<double>f2)->bool{return (f1[1]<f2[1]);});

        int lambda_i = 0;
        int mu_i = 0;
        // populate with the best individuals of both populations
        for (int i = 0 ; i < parent_population.mu; i++) {
            if(fitness_lambda[lambda_i][1] <= fitness_mu[mu_i][1]) {
                parent_population.individuals[i] = lambda_population.individuals[fitness_lambda[lambda_i][0]];
                lambda_i++;
            } else {
                parent_population.individuals[i] = parent_population.individuals[fitness_mu[mu_i][0]];
                mu_i++;
            }
        }

    }


}





















