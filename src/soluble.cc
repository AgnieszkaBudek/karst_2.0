/*
 * soluble.cc
 *
 *  Created on: May 18, 2020
 *      Author: agnieszka
 */

#include "soluble.h"

#include "solid.h"
#include "pore.h"
#include "node.h"
#include "network.h"

#include <math.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>
#include <string.h>


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

	double f     = p->local_Da_eff(S);      //effective reaction rate (taking into account both reaction and transversal diffusion)
	double q_tmp = fabs(p->q);

	return q_tmp*exp(-f);

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


	double dd       = S->R->default_diameter_change(p0,S);
	double q_tmp    = fabs(p0->q);
	double c_tmp_in = p0->calculate_inlet_c(0);


	double x = 0;   //value to be returned

	if(p0->d + dd < S->d_min)   { //if there is no space for full precipitation
		double dd_plus = S->R->m[0]->dd(p0,S);
		double dd_minus=(p0->d/S->d0 + dd_plus - S->d_min/S->d0);
		x = c_tmp_in*q_tmp*(1-exp(-f1)) - (M_PI*(p0->d)*(dd_minus*S->d0)/2*p0->l)/(S->gamma*S->dt/S->dt_unit);
		if(x<0) {
			if (S->if_verbose) cerr<< "WARNING: STH wrong in outlet_c_c_1 with calculating x!!!"<<endl;
			x=0;}
		}
	else{	//if there is no problem with a space for precipitation
		if (f1!=f2) x = c_tmp_in*q_tmp*    (exp(-f1) - exp(-f2))*(f1)/(f2-f1);
		else        x = c_tmp_in*q_tmp*    f2*exp(-f2);
		}

	return x;

}

double Soluble_C_Agnieszkas_model::c1(Pore *p0, Network *S){


	if(!is_reacting(p0,S)) return 1;

	double f2       = p0->local_Da_eff_2 (S);
	double dd       = S->R->default_diameter_change(p0,S);
	double q_tmp    = fabs(p0->q);

	//Checking if there is enough space for full dissolution
	if(p0->d + dd < S->d_min){
		double delta_minus = S->R->m[1]->dd(p0,S);
		if(p0->is_V_left(0)) return q_tmp;            //if there is no space for full precipitation I use ugly but working formula form outlet_c_c_2
		else                 return q_tmp*(1 - (p0->d-S->d_min)/fabs(delta_minus)*(1-exp(-f2)));
	}
	else 		return  q_tmp*exp(-f2);	  //if there is enough space for precipitation

}

