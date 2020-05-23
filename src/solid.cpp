/*
 * solid.cpp
 *
 *  Created on: May 18, 2020
 *      Author: agnieszka
 */

#include "solid.h"
#include "pore.h"
#include "node.h"
#include "network.h"

Solid::Solid() {
	// TODO Auto-generated constructor stub

}

Solid::~Solid() {
	// TODO Auto-generated destructor stub
}

bool Solid_E_Agnieszkas_model::is_reacting(Pore *p, Network *S){
	if(p->d==0 || p->q==0)                   return false;
	if(p->l==S->l_min)                       return false;   //no reaction in tiny grain
	if(p->d<=S->d_min && (!p->is_V_left(0))) return false;

	return true;
}

double Solid_E_Agnieszkas_model::dd(Pore *p, Network *S){

	if (!is_reacting(p,S)) return 0;

	//dissolution parameters
	double f1      = p->local_Da_eff(S);
	double g       = p->local_G(S);
	double c0      = p->calculate_inlet_c(0);


	//precipitation parameters
	double f2       = p->local_Da_eff_2(S);
	double c0_c     = p->calculate_inlet_c(1);
	double dd_minus = 0; 		//diameter change


	//finding precipitation contribution
	if      (f2==0)         			 dd_minus = 0;
	else if (f1==f2 && p->is_V_left(0))  dd_minus = S->gamma*S->dt/(1+g)/f1*((c0_c + c0)*(1-exp(-f1)) -c0*exp(-f1)*f1);
	else if (!p->is_V_left(0)) 		     dd_minus = S->gamma*S->dt/(1+g)/f1*  c0_c*      (1-exp(-f2));
	else                    			 dd_minus = S->gamma*S->dt/(1+g)/f1*(\
												   c0  * (f1*(1-exp(-f2)) - f2*(1-exp(-f1)))/(f1-f2)+\
												   c0_c*  (1-exp(-f2)) );

	return dd_minus;
}

bool Solid_A_Agnieszkas_model::is_reacting(Pore *p, Network *S){

	if(p->d==0||p->q==0)                      return false;
	if(p->l<=S->l_min)                        return false;
	if(S->if_track_grains &&!p->is_V_left(0)) return false;

	return true;
}

double Solid_A_Agnieszkas_model::dd(Pore *p, Network *S){

	if(!is_reacting(p,S)) return 0;

	//dissolution parameters
	double f1      = p->local_Da_eff(S);
	double g       = p->local_G(S);
	double c0;
	c0 = p->calculate_inlet_c(0);

	double dd_plus = 0; 		//diameter change


	//finding dissolution contribution
	if      (f1==0)      dd_plus = 0;
	else if (S->G1 >=0)  dd_plus = S->dt*c0*(1-exp(-f1))/(1+g)/f1;
	else        	     dd_plus = S->dt*c0*(1-exp(-f1))/f1/p->d;


	return dd_plus;

}
