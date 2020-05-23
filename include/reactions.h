/*
 * reactions.h
 *
 *  Created on: May 17, 2020
 *      Author: agnieszka
 */

#ifndef INCLUDE_REACTIONS_H_
#define INCLUDE_REACTIONS_H_

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>

#include "soluble.h"
#include "solid.h"

class Soluble;
class Solid;
class Pore;

class Reactions {

public:

	int bw;    ///< number of soluble species
	int bm;    ///< number of solid species
	int br;    ///< number of independent reactions

	Soluble** w;    	///< table with soluble materials
	Solid  ** m;   		///< table with solid materials

	double* k;     		///< table with reactions rate (not used now)
	double* D;     		///< table with diffusion coefficient (not used now)
	double* DD;    		///< table with transversal diffusion coefficient (not used now)
	double* gamma; 		///< table with capacity number for all solid materials
	double* C;			///< inlet concentration of all soluble materials
	double* V;		    ///< table with total amount of volume for all species


public:
	// Different models that are implemented:
	void Set_reactions_to_pure_dissolution();  						///< Model with pure dissolution due to linear reaction with acid injected to the sysytem
	void Set_reactions_to_dissolution_and_precipitation_A();        ///< Model with both dissolution and precipitation (my model form PhD studies)
	void Set_reactions_to_dissolution_and_precipitation_F1();       ///< First approach to simulate Florians model

	Reactions(string model_name);
	virtual ~Reactions();



};

#endif /* INCLUDE_REACTIONS_H_ */
