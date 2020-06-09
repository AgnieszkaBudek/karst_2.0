/*
 * soluble.h
 *
 *  Created on: May 18, 2020
 *      Author: agnieszka
 */

#ifndef INCLUDE_SOLUBLE_H_
#define INCLUDE_SOLUBLE_H_

#include "reactions.h"

class Pore;
class Network;
class Singla_Reaction;

class Soluble {

public:

	double D;           ///< radial diffusion coefficient
	double DD;          ///< transversal diffusion coefficient
	double C0;          ///< concentration at the inlet

	int br;              ///< list of reactions it takes part //not used
	Single_Reaction **r; ///< list of reactions it takes part //not used

	virtual double  c1(Pore *p, Network *S)=0;           //c(i+1) = c1*c(i) + c0 ; coefficient in linear eq for concentration
	virtual double  c0(Pore *p, Network *S)=0;           //c(i+1) = c1*c(i) + c0 ; constant in linear eq. for concentration

	virtual bool   is_reacting(Pore *p, Network *S)=0;

public:
	Soluble();
	Soluble(int bbr, Single_Reaction** rr);
	virtual ~Soluble();
};



class Soluble_B_Agnieszkas_model : public Soluble {

	public:
	virtual double c1(Pore *p, Network *S);
	virtual double c0(Pore *p, Network *S);

	virtual bool   is_reacting(Pore *p, Network *S);

};



class Soluble_C_Agnieszkas_model : public Soluble {

	public:
	virtual double c1(Pore *p, Network *S);
	virtual double c0(Pore *p, Network *S);

	virtual bool   is_reacting(Pore *p, Network *S);

};

#endif /* INCLUDE_SOLUBLE_H_ */
