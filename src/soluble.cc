/*
 * soluble.cc
 *
 *  Created on: May 18, 2020
 *      Author: agnieszka
 */

#include "../include/soluble.h"

#include "solid.h"
#include "pore.h"
#include "node.h"
#include "network.h"
#include <math.h>

Soluble::Soluble() {
	// TODO Auto-generated constructor stub

}

Soluble::~Soluble() {
	// TODO Auto-generated destructor stub
}

bool Soluble_B_Agnieszkas_model::is_reacting(Pore *p, Network *S){

	if(p->d==0||p->q==0)                      return false;
	if(p->l<=S->l_min)                        return false;
	if(S->if_track_grains &&!p->is_V_left(0)) return false;

	return true;
}

double Soluble_B_Agnieszkas_model::c0(Pore *p, Network *S){
	return 0;
}

double Soluble_B_Agnieszkas_model::c1(Pore *p, Network *S){

	if(!is_reacting(p,S)) return 1;

	double f = p->local_Da_eff(S);      //effective reaction rate (taking into account both reaction and transversal diffusion)

	return exp(-f);

}


bool Soluble_C_Agnieszkas_model::is_reacting(Pore *p, Network *S){

	if(p->d==0 || p->q==0)                   return false;
	if(p->l==S->l_min)                       return false;   //no reaction in tiny grain
	if(p->d<=S->d_min && (!p->is_V_left(0))) return false;

	return true;
}


double Soluble_C_Agnieszkas_model::c0(Pore *p0, Network *S){

	if(!is_reacting(p0,S))    return 0;   //no contribution if there is no reaction
	if(!(p0->is_V_left(0)))   return 0;   //no contribution if there is no A material

	double f1 = p0->local_Da_eff   (S);
	double f2 = p0->local_Da_eff_2 (S);

	if (f1!=f2) return   (exp(-f1) - exp(-f2))*(f1)/(f2-f1);
	else        return    f2*exp(-f2);

}

double Soluble_C_Agnieszkas_model::c1(Pore *p0, Network *S){

	if(!is_reacting(p0,S)) return 1;

	double f2       = p0->local_Da_eff_2 (S);

	return  exp(-f2);
}

