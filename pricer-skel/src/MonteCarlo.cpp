/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "MonteCarlo.hpp";
#include "BlackScholesModel.hpp";
#include "Option.hpp";
#include <math.h>;

    /**
     * Calcule le prix de l'option à la date 0
     *
     * @param[out] prix valeur de l'estimateur Monte Carlo
     * @param[out] ic largeur de l'intervalle de confiance
     */
    void price(double &prix, double &ic) {
        PnlMat * path = pnl_mat_create(this->opt_->size_, this->opt_->nbTimeSteps_);
        PnlVect * payoff_montercarlo_vect = pnl_vect_create(this->nbSamples_);
        
        for(int i=0; i< this->nbSamples_; i++) {
           this->mod_->asset(path,this->opt_->T_,this->opt_->nbTimeSteps_,this->rng_);
           pnl_vect_set(payoff_montercarlo_vect, i, this->opt_->payoff(path));
        }
        
        double average = 0;
        double square_average = 0;
        
        for(int i=0; i<this->nbSamples_; i++){
            double payoff_i = pnl_vect_get(payoff_montercarlo_vect,i);
            average += payoff_i;
            square_average += pow(payoff_i); 
        }
        average /= this->nbSamples_;
        square_average /= this->nb_Samples_;
        
        double varianceMontecarlo =  exp(-2 * this->opt_->T_ * this->mod_->r_) *
            (square_average - pow(average,2));
        
        
        prix = average * exp(this->opt_->T_ * this->mod_->r_);
        ic = 1.96 * 2 * varianceMontecarlo / sqrt(this->nb_Samples_);
    }


